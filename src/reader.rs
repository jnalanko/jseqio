use ex::fs::File; // File streams that include the filename in the error messages
use std::io::{BufRead, BufReader};
use std::path::Path;
use flate2::read::MultiGzDecoder;
use crate::{FileType};
use crate::record::RefRecord;

// Takes a BufRead because we need read_until.
pub struct FastXReader<R: std::io::BufRead>{
    pub filetype: FileType,
    pub filename: Option<String>, // Only used for error messages. If None, the file is unknown or there is no file, like when reading from stdin.
    pub input: R,
    pub seq_buf: Vec<u8>,
    pub head_buf: Vec<u8>,
    pub qual_buf: Vec<u8>,
    pub plus_buf: Vec<u8>, // For the fastq plus-line
    pub fasta_temp_buf: Vec<u8>, // Stores the fasta header read in the previous iteration
}

// Trait for a stream returning SeqRecord objects, used in DynamicFastXReader to abstract over
// The input stream type.
trait SeqRecordProducer {
    fn read_next(&mut self) -> Result<Option<RefRecord>, Box<dyn std::error::Error>>;

    // Since we want to call this for trait objects where we don't know the size of the struct,
    // We need to take self in a Box.
    fn into_db_boxed(self: Box<Self>) -> Result<crate::seq_db::SeqDB, Box<dyn std::error::Error>>;

    fn filetype(&self)-> FileType; 

    // For error messages
    fn set_filepath(&mut self, filepath: &Path);
}

#[derive(Debug)]
pub struct ParseError{
    pub message: String,
    pub filename: Option<String>,
    pub filetype: Option<FileType>,
}

impl std::error::Error for ParseError{}

impl std::fmt::Display for ParseError{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl<R: std::io::BufRead> FastXReader<R>{

    fn build_parse_error(&self, message: &str) -> Box<ParseError>{
        Box::new(
            ParseError{
                message: message.to_owned(), 
                filename: self.filename.clone(), 
                filetype: Some(self.filetype)
            }
        )
    }

    pub fn read_next(&mut self) -> Result<Option<RefRecord>, Box<dyn std::error::Error>> {
        if matches!(self.filetype, FileType::FASTQ){
            // FASTQ format

            self.seq_buf.clear();
            self.head_buf.clear();
            self.qual_buf.clear();
            self.plus_buf.clear();

            // Read header line
            let bytes_read = self.input.read_until(b'\n', &mut self.head_buf)?;
            if bytes_read == 0 {return Ok(None)} // End of stream
            if self.head_buf[0] != b'@'{
                return Err(self.build_parse_error("FASTQ header line does not start with @"));
            }

            // Read sequence line
            let bytes_read = self.input.read_until(b'\n', &mut self.seq_buf)?;
            if bytes_read == 0 {
                return Err(self.build_parse_error("FASTQ sequence line missing.")); // File can't end here
            }
            
            // read +-line
            let bytes_read = self.input.read_until(b'\n', &mut self.plus_buf)?;
            if bytes_read == 0 {
                return Err(self.build_parse_error("FASTQ + line missing.")); // File can't end here
            }

            // read qual-line
            let bytes_read = self.input.read_until(b'\n', &mut self.qual_buf)?;
            if bytes_read == 0 { // File can't end here
                return Err(self.build_parse_error("FASTQ quality line missing.")); // File can't end here
            } else if bytes_read != self.seq_buf.len(){
                let msg = format!("FASTQ quality line has different length than sequence line ({} vs {})", bytes_read, self.seq_buf.len());
                return Err(self.build_parse_error(&msg));
            }

            Ok(Some(RefRecord{head: self.head_buf.as_slice().strip_prefix(b"@").unwrap().strip_suffix(b"\n").unwrap(), 
                                seq: self.seq_buf.as_slice().strip_suffix(b"\n").unwrap(),
                                qual: Some(self.qual_buf.as_slice().strip_suffix(b"\n").unwrap())}))
        }
        else{
            // FASTA format
            self.seq_buf.clear();
            self.head_buf.clear();

            // Read header line
            if self.fasta_temp_buf.is_empty() {
                // This is the first record -> read header from input
                let bytes_read = self.input.read_until(b'\n', &mut self.head_buf)?;
                if bytes_read == 0 {return Ok(None)} // End of stream
            } else{
                // Take stashed header from previous iteration
                self.head_buf.append(&mut self.fasta_temp_buf); // Also clears the temp buf
            }

            // Read sequence line
            loop{
                let bytes_read = self.input.read_until(b'\n', &mut self.fasta_temp_buf)?;
                if bytes_read == 0 {
                    // No more bytes left to read
                    if self.seq_buf.is_empty(){
                        // Stream ends with an empty sequence
                        return Err(self.build_parse_error("Empty sequence in FASTA file"));
                    }
                    break; // Ok, last record of the file
                }

                // Check if we read the header of the next read
                let start = self.fasta_temp_buf.len() as isize - bytes_read as isize;
                if self.fasta_temp_buf[start as usize] == b'>'{
                    // Found a header. Leave it to the buffer for the next iteration.
                    break;
                } else{
                    // Found more sequence -> Append to self.seq_buf
                    self.seq_buf.append(&mut self.fasta_temp_buf); // Also clears the temp buf
                    self.seq_buf.pop(); // Trim newline (TODO: what if there is none?)
                }
            }
            Ok(Some(RefRecord{head: self.head_buf.as_slice().strip_prefix(b">").unwrap().strip_suffix(b"\n").unwrap(), 
                                seq: self.seq_buf.as_slice(), // Newlines are already trimmed before
                                qual: None}))
        }
    }

    pub fn new_with_format(input: R, filetype: FileType) -> Self{
        FastXReader{filetype,
                    input,
                    filename: None,
                    seq_buf: Vec::<u8>::new(),
                    head_buf: Vec::<u8>::new(),
                    qual_buf: Vec::<u8>::new(),
                    plus_buf: Vec::<u8>::new(),
                    fasta_temp_buf: Vec::<u8>::new(),}
    }

    // Detect whether it's fasta or FASTQ from the first byte.
    pub fn new(mut input: R) -> Result<Self, Box<dyn std::error::Error>>{
        let bytes = input.fill_buf()?;

        // Empty file is arbitrarily considered as FASTA
        let mut filetype = FileType::FASTA;

        if !bytes.is_empty(){
            filetype = match bytes[0]{
                b'>' => FileType::FASTA,
                b'@' => FileType::FASTQ,
                _ => return Err(
                        Box::new(ParseError{message: "Error: File does not start with '>' or '@'".to_owned(), 
                        filename: None, 
                        filetype: None}))
            } 
        }

        Ok(FastXReader{filetype,
            input,
            filename: None,
            seq_buf: Vec::<u8>::new(),
            head_buf: Vec::<u8>::new(),
            qual_buf: Vec::<u8>::new(),
            plus_buf: Vec::<u8>::new(),
            fasta_temp_buf: Vec::<u8>::new()})

    }

    pub fn into_db(mut self) -> Result<crate::seq_db::SeqDB, Box<dyn std::error::Error>>{
        let mut headbuf: Vec<u8> = Vec::new();
        let mut seqbuf: Vec<u8> = Vec::new();
        let mut qualbuf: Option<Vec<u8>> = match self.filetype{
            FileType::FASTQ => Some(Vec::new()),
            FileType::FASTA => None,
        };

        let mut head_starts: Vec<usize> = vec![0];
        let mut seq_starts: Vec<usize> = vec![0];
        let mut qual_starts: Option<Vec<usize>> = match self.filetype{
            FileType::FASTQ => Some(vec![0]),
            FileType::FASTA => None,
        };

        while let Some(rec) = self.read_next()?{
            headbuf.extend_from_slice(rec.head);
            seqbuf.extend_from_slice(rec.seq);
            if let Some(qual) = rec.qual{
                qualbuf.as_mut().expect("Error: found a fastq record in a fasta stream.").extend_from_slice(qual);
                qual_starts.as_mut().unwrap().push(qualbuf.as_ref().unwrap().len());
            }
            head_starts.push(headbuf.len());
            seq_starts.push(seqbuf.len());
        }
        Ok(crate::seq_db::SeqDB{headbuf, seqbuf, qualbuf, head_starts, seq_starts, qual_starts})
    }



}


pub struct DynamicFastXReader {
    stream: Box<dyn SeqRecordProducer>,
}

// A class that contains a dynamic trait object for different
// types of input streams.
impl DynamicFastXReader {
 
    // New from file
    pub fn from_file<P: AsRef<std::path::Path>>(filepath: &P) -> Result<Self, Box<dyn std::error::Error>> {
        let input = File::open(filepath).unwrap();
        let mut reader = Self::new(BufReader::new(input))?;
        reader.stream.set_filepath(filepath.as_ref());
        Ok(reader)
    }

    // New from stdin
    pub fn from_stdin() -> Result<Self, Box<dyn std::error::Error>> {
        let input = std::io::stdin();
        let reader = Self::new(BufReader::new(input))?;
        Ok(reader)
    }

    // New from stream, with automatic gzip detection
    pub fn new<R: std::io::BufRead + 'static>(mut input: R) -> Result<Self, Box<dyn std::error::Error>>{
        let bytes = input.fill_buf()?;
        let mut gzipped = false;
        match bytes.len(){
            0 => (), // Empty file
            1 => return Err(Box::new(ParseError{message: "Corrupt FASTA/FASTQ file: only one byte found.".to_owned(), 
                        filename: None, 
                        filetype: None})),
            _ => { // Two or more bytes available. Check if the first two are a valid gzip header.
                if bytes[0] == 0x1f && bytes[1] == 0x8b{ 
                    gzipped = true;
                }
            }
        }

        match gzipped{
            true => {
                let gzdecoder = MultiGzDecoder::<R>::new(input);

                // We wrap this in BufReader because the FastX parser requires buffered reading
                let gzbufdecoder = BufReader::<MultiGzDecoder::<R>>::new(gzdecoder);
                Self::from_raw_stream(gzbufdecoder)
            },
            false => Self::from_raw_stream(input)
        }
    }

    // Creates a reader from a raw stream of uncompressed data (no gzip detection). Used by other constructors.
    // Need to constrain + 'static because boxed trait objects always need to have a static lifetime.
    fn from_raw_stream<R: std::io::BufRead + 'static>(r: R) -> Result<Self, Box<dyn std::error::Error>>{
        let reader = FastXReader::<R>::new(r)?;
        Ok(DynamicFastXReader {stream: Box::new(reader)})
    }

    // Returns None if no more records
    pub fn read_next(&mut self) -> Result<Option<RefRecord>, Box<dyn std::error::Error>>{
        self.stream.read_next()
    }

    pub fn filetype(&self)-> FileType{
        self.stream.filetype()
    }

    pub fn into_db(self) -> Result<crate::seq_db::SeqDB, Box<dyn std::error::Error>>{
        self.stream.into_db_boxed()
    }

}

// Implement common SeqStream trait for all
// FastXReaders over the generic parameter R.
impl<R: BufRead> SeqRecordProducer for FastXReader<R>{

    fn read_next(&mut self) -> Result<Option<RefRecord>, Box<dyn std::error::Error>>{
        self.read_next()
    }

    fn filetype(&self)-> FileType{
        self.filetype
    }

    fn into_db_boxed(self: Box<Self>) -> Result<crate::seq_db::SeqDB, Box<dyn std::error::Error>>{
        self.into_db()
    }

    fn set_filepath(&mut self, filepath: &Path){
        self.filename = Some(filepath.as_os_str().to_str().unwrap().to_owned());
    }

}

