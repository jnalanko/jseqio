use ex::fs::File; // File streams that include the filename in the error messages
use std::io::{BufRead, BufReader};
use std::path::Path;
use flate2::read::MultiGzDecoder;
use crate::seq_db::SeqDB;
use crate::{FileType};
use crate::record::{MutRefRecord, OwnedRecord, RefRecord};

// Takes a BufRead because we need read_until.
pub struct StaticFastXReader<R: std::io::BufRead>{
    pub filetype: FileType,
    pub filename: Option<String>, // Only used for error messages. If None, the file is unknown or there is no file, like when reading from stdin.
    pub input: R,
    pub seq_buf: Vec<u8>,
    pub head_buf: Vec<u8>,
    pub qual_buf: Vec<u8>,
    pub plus_buf: Vec<u8>, // For the fastq plus-line
    pub fasta_temp_buf: Vec<u8>, // Stores the fasta header read in the previous iteration
}

pub trait SeqStream{
    fn read_next(&mut self) -> Result<Option<RefRecord>, Box<dyn std::error::Error>>;
}

// Trait for a stream returning RefRecord objects, used in DynamicFastXReader to abstract over
// The input stream type.
trait JSeqIOReaderInterface{
    fn read_next(&mut self) -> Result<Option<RefRecord>, Box<dyn std::error::Error>>;
    fn read_next_mut(&mut self) -> Result<Option<MutRefRecord>, Box<dyn std::error::Error>>;

    // Since we want to call this for trait objects where we don't know the size of the struct,
    // We need to take self in a Box.
    fn into_db_boxed(self: Box<Self>) -> Result<crate::seq_db::SeqDB, Box<dyn std::error::Error>>;
    fn into_db_with_revcomp_boxed(self: Box<Self>) -> Result<(crate::seq_db::SeqDB, crate::seq_db::SeqDB), Box<dyn std::error::Error>>;

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

impl<R: std::io::BufRead> StaticFastXReader<R>{

    fn build_parse_error(&self, message: &str) -> Box<ParseError>{
        Box::new(
            ParseError{
                message: message.to_owned(), 
                filename: self.filename.clone(), 
                filetype: Some(self.filetype)
            }
        )
    }

    fn read_fasta_record(&mut self) -> Result<Option<MutRefRecord>, Box<dyn std::error::Error>>{
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

        // Make sure all letters are in same case
        for c in self.seq_buf.iter_mut() {
            c.make_ascii_uppercase();
        }

        // Trim '>' and '\n' from the header 
        assert!(self.head_buf[0] == b'>');
        assert!(*self.head_buf.last().unwrap() == b'\n');
        let head_buf_len = self.head_buf.len();
        let head = &mut self.head_buf.as_mut_slice()[1..head_buf_len-1]; // Remove '>' and '\n'

        Ok(Some(MutRefRecord{head: head, 
                            seq: self.seq_buf.as_mut_slice(), // Newlines are already trimmed before
                            qual: None}))
    }

    fn read_fastq_record(&mut self) -> Result<Option<MutRefRecord>, Box<dyn std::error::Error>>{

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
        assert!(self.head_buf.last().unwrap() == &b'\n');
        self.head_buf.pop(); // Trim the '\n'

        // Read sequence line
        let bytes_read = self.input.read_until(b'\n', &mut self.seq_buf)?;
        if bytes_read == 0 {
            return Err(self.build_parse_error("FASTQ sequence line missing.")); // File can't end here
        }
        assert!(self.seq_buf.last().unwrap() == &b'\n');
        self.seq_buf.pop(); // Trim the '\n'
        
        // read +-line
        let bytes_read = self.input.read_until(b'\n', &mut self.plus_buf)?;
        if bytes_read == 0 {
            return Err(self.build_parse_error("FASTQ + line missing.")); // File can't end here
        }

        // read qual-line
        let bytes_read = self.input.read_until(b'\n', &mut self.qual_buf)?;
        if bytes_read == 0 { // File can't end here
            return Err(self.build_parse_error("FASTQ quality line missing.")); // File can't end here
        } 
        assert!(self.qual_buf.last().unwrap() == &b'\n');
        self.qual_buf.pop(); // Trim the '\n'
        if self.qual_buf.len() != self.seq_buf.len() {
            let msg = format!("FASTQ quality line has different length than sequence line ({} vs {})", self.qual_buf.len(), self.seq_buf.len());
            return Err(self.build_parse_error(&msg));
        }

        // Make sure all letters are in same case
        for c in self.seq_buf.iter_mut() {
            c.make_ascii_uppercase();
        }

        assert!(self.head_buf[0] == b'@');
        Ok(Some(MutRefRecord{head: &mut self.head_buf[1..], // Remove '@'
                            seq: &mut self.seq_buf,
                            qual: Some(&mut self.qual_buf)}))
    }

    // Read one record from the input.
    // This is not named just next() because it's not a Rust iterator because it streams the input.
    pub fn read_next(&mut self) -> Result<Option<RefRecord>, Box<dyn std::error::Error>> {
        match self.filetype{
            FileType::FASTA => self.read_fasta_record().map(|opt| opt.map(|rec| rec.into_shared_ref())),
            FileType::FASTQ => self.read_fastq_record().map(|opt| opt.map(|rec| rec.into_shared_ref())),
        }
    }

    // Read one record from the input.
    // This is not named just next() because it's not a Rust iterator because it streams the input.
    pub fn read_next_mut(&mut self) -> Result<Option<MutRefRecord>, Box<dyn std::error::Error>> {
        match self.filetype{
            FileType::FASTA => self.read_fasta_record(),
            FileType::FASTQ => self.read_fastq_record(),
        }
    }

    // New with known format
    fn new_with_format(input: R, filetype: FileType) -> Self{
        StaticFastXReader{filetype,
                    input,
                    filename: None,
                    seq_buf: Vec::<u8>::new(),
                    head_buf: Vec::<u8>::new(),
                    qual_buf: Vec::<u8>::new(),
                    plus_buf: Vec::<u8>::new(),
                    fasta_temp_buf: Vec::<u8>::new(),}
    }

    // Detect whether it's fasta or FASTQ based on the first byte.
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

        Ok(StaticFastXReader::new_with_format(input, filetype))

    }

    pub fn into_db_with_revcomp(mut self) -> Result<(SeqDB, SeqDB), Box<dyn std::error::Error>>{

        // Reusable record for storing the reverse complement
        let mut rc_record = 
        OwnedRecord{
            head: Vec::new(), 
            seq: Vec::new(),
            qual: match self.filetype{
                FileType::FASTA => None,
                FileType::FASTQ => Some(Vec::new()),
            }
        };

        let mut fw_db = SeqDB::new();
        let mut rc_db = SeqDB::new();

        while let Some(rec) = self.read_next()?{
            fw_db.push_record(rec);

            // Copy fw record data to rc record
            rc_record.head.clear();
            rc_record.head.extend_from_slice(rec.head);
            rc_record.seq.clear();
            rc_record.seq.extend_from_slice(rec.seq);
            if let Some(qual) = &mut rc_record.qual{
                qual.clear();
                qual.extend_from_slice(rec.qual.unwrap());
            }

            // Reverse complement and push to database
            rc_record.reverse_complement();
            rc_db.push_record(rc_record.as_ref_record());
        }
        fw_db.shrink_to_fit();
        rc_db.shrink_to_fit();
        Ok((fw_db, rc_db))
    }

    pub fn into_db(mut self) -> Result<crate::seq_db::SeqDB, Box<dyn std::error::Error>>{
        let mut db = SeqDB::new();

        while let Some(rec) = self.read_next()?{
            db.push_record(rec);
        }
        db.shrink_to_fit();
        Ok(db)
    }



}


pub struct DynamicFastXReader {
    stream: Box<dyn JSeqIOReaderInterface + Send>,
    compression_type: crate::CompressionType,
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
    pub fn new<R: std::io::BufRead + 'static + Send>(mut input: R) -> Result<Self, Box<dyn std::error::Error>>{
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
                Self::from_raw_stream(gzbufdecoder, crate::CompressionType::Gzip)
            },
            false => Self::from_raw_stream(input, crate::CompressionType::None)
        }
    }

    // Creates a reader from a raw stream of uncompressed data (no gzip detection). Used by other constructors.
    // Need to constrain + 'static because boxed trait objects always need to have a static lifetime.
    fn from_raw_stream<R: std::io::BufRead + 'static + Send>(r: R, compression_type: crate::CompressionType) -> Result<Self, Box<dyn std::error::Error>>{
        let reader = StaticFastXReader::<R>::new(r)?;
        Ok(DynamicFastXReader {stream: Box::new(reader), compression_type})
    }

    pub fn into_db(self) -> Result<crate::seq_db::SeqDB, Box<dyn std::error::Error>>{
        self.stream.into_db_boxed()
    }

    pub fn into_db_with_revcomp(self) -> Result<(crate::seq_db::SeqDB, crate::seq_db::SeqDB), Box<dyn std::error::Error>>{
        self.stream.into_db_with_revcomp_boxed()
    }

    pub fn compression_type(&self) -> crate::CompressionType{
        self.compression_type
    }

    pub fn read_next(&mut self) -> Result<Option<RefRecord>, Box<dyn std::error::Error>>{
        self.stream.read_next()
    }

    pub fn read_next_mut(&mut self) -> Result<Option<MutRefRecord>, Box<dyn std::error::Error>>{
        self.stream.read_next_mut()
    }

    pub fn filetype(&self)-> FileType{
        self.stream.filetype()
    }
    
    // For error messages
    pub fn set_filepath(&mut self, filepath: &Path){
        self.stream.set_filepath(filepath);
    }
}


// Implement common SeqRecordProducer trait for all
// StaticFastXReaders over the generic parameter R.
impl<R: BufRead> JSeqIOReaderInterface for StaticFastXReader<R>{

    fn read_next(&mut self) -> Result<Option<RefRecord>, Box<dyn std::error::Error>>{
        self.read_next()
    }

    fn read_next_mut(&mut self) -> Result<Option<MutRefRecord>, Box<dyn std::error::Error>> {
        self.read_next_mut()
    }

    fn filetype(&self)-> FileType{
        self.filetype
    }

    fn into_db_boxed(self: Box<Self>) -> Result<crate::seq_db::SeqDB, Box<dyn std::error::Error>>{
        self.into_db()
    }

    fn into_db_with_revcomp_boxed(self: Box<Self>) -> Result<(crate::seq_db::SeqDB, crate::seq_db::SeqDB), Box<dyn std::error::Error>>{
        self.into_db_with_revcomp()
    }    

    // For error messages
    fn set_filepath(&mut self, filepath: &Path){
        self.filename = Some(filepath.as_os_str().to_str().unwrap().to_owned());
    }

}

impl SeqStream for DynamicFastXReader {
    fn read_next(&mut self) -> Result<Option<RefRecord>, Box<dyn std::error::Error>> {
        DynamicFastXReader::read_next(self)
    }
}

impl<R: BufRead> SeqStream for StaticFastXReader<R> {
    fn read_next(&mut self) -> Result<Option<RefRecord>, Box<dyn std::error::Error>> {
        StaticFastXReader::read_next(self)
    }
}

// Turns a SeqStream into another SeqStream that also streams the reverse
// complements of the original stream. See read_next() for more information.
pub struct SeqStreamWithRevComp<S: SeqStream> {
    inner: S,
    rec: OwnedRecord,
    parity: bool, // Every other sequence we return is a reverse complement of the previous
}

impl<S: SeqStream> SeqStreamWithRevComp<S> {
    pub fn new(inner: S) -> Self{
        Self{
            inner,
            rec: OwnedRecord{seq: Vec::new(), head: Vec::new(), qual: None},
            parity: false,
        }
    }
}

impl<S: SeqStream> SeqStream for SeqStreamWithRevComp<S> {
    // If the original sequence stream is x1,x2,x3..., returns sequences in the order
    // x1, rc(x1), x2, rc(x2), x3, rc(x3)...
    fn read_next(&mut self) -> Result<Option<RefRecord>, Box<dyn std::error::Error>> {
        self.parity = !self.parity;

        if self.parity {
            let new = match self.inner.read_next()? {
                None => return Ok(None), // End of stream
                Some(r) => r
            };

            // Copy the record to self storage for reverse complementation later
            self.rec.head.clear();
            self.rec.head.extend(new.head);
            self.rec.seq.clear();
            self.rec.seq.extend(new.seq);
            if let Some(q) = new.qual {
                if self.rec.qual.is_none() {
                    self.rec.qual = Some(Vec::<u8>::new());
                }
                self.rec.qual.as_mut().unwrap().clear();
                self.rec.qual.as_mut().unwrap().extend(q);
            }

            return Ok(Some(self.rec.as_ref_record()));

        } else {
            crate::reverse_complement_in_place(&mut self.rec.seq);
            return Ok(Some(self.rec.as_ref_record()));
        }
    }
}
