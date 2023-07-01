use std::fs::File;
use std::io::{BufRead, BufReader};
use flate2::read::MultiGzDecoder;
use crate::{FileType, figure_out_file_format};
use crate::record::RefRecord;

// Takes a BufRead because we need read_until.
pub struct FastXReader<R: std::io::BufRead>{
    pub filetype: FileType,
    pub input: R,
    pub seq_buf: Vec<u8>,
    pub head_buf: Vec<u8>,
    pub qual_buf: Vec<u8>,
    pub plus_buf: Vec<u8>, // For the fastq plus-line
    pub fasta_temp_buf: Vec<u8>, // Stores the fasta header read in the previous iteration
}

impl<R: std::io::BufRead> FastXReader<R>{
    pub fn read_next(&mut self) -> Option<RefRecord>{
        if matches!(self.filetype, FileType::FASTQ){
            // FASTQ format

            self.seq_buf.clear();
            self.head_buf.clear();
            self.qual_buf.clear();
            self.plus_buf.clear();

            // Read header line
            let bytes_read = self.input.read_until(b'\n', &mut self.head_buf);
            if bytes_read.expect("I/O error.") == 0 {return None} // End of stream

            // Read sequence line
            let bytes_read = self.input.read_until(b'\n', &mut self.seq_buf);
            if bytes_read.expect("I/O error.") == 0 {
                panic!("FASTQ sequence line missing."); // File can't end here
            }
            
            // read +-line
            let bytes_read = self.input.read_until(b'\n', &mut self.plus_buf);
            if bytes_read.expect("I/O error.") == 0 {
                panic!("FASTQ + line missing."); // File can't end here
            }

            // read qual-line
            let bytes_read = self.input.read_until(b'\n', &mut self.qual_buf);
            let bytes_read = bytes_read.expect("I/O error.");
            if bytes_read == 0{ // File can't end here
                panic!("FASTQ quality line missing."); 
            } else if bytes_read != self.seq_buf.len(){
                panic!("FASTQ quality line has different length than sequence line ({} vs {})", bytes_read, self.seq_buf.len())
            }

            return Some(RefRecord{head: self.head_buf.as_slice().strip_prefix(b"@").unwrap().strip_suffix(b"\n").unwrap(), 
                                seq: self.seq_buf.as_slice().strip_suffix(b"\n").unwrap(),
                                qual: Some(self.qual_buf.as_slice().strip_suffix(b"\n").unwrap())})
        }
        else{
            // FASTA format
            self.seq_buf.clear();
            self.head_buf.clear();

            // Read header line
            if self.fasta_temp_buf.is_empty() {
                // This is the first record -> read header from input
                let bytes_read = self.input.read_until(b'\n', &mut self.head_buf);
                if bytes_read.expect("I/O error.") == 0 {return None} // End of stream
            } else{
                // Take stashed header from previous iteration
                self.head_buf.append(&mut self.fasta_temp_buf); // Also clears the temp buf
            }

            // Read sequence line
            loop{
                let bytes_read = self.input.read_until(b'\n', &mut self.fasta_temp_buf);
                match bytes_read{
                    Err(e) => panic!("{}",e), // File can't end here
                    Ok(bytes_read) => {
                        if bytes_read == 0{
                            // No more bytes left to read
                            if self.seq_buf.is_empty(){
                                // Stream ends with an empty sequence
                                panic!("Empty sequence in FASTA file");
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
                }
            }

            return Some(RefRecord{head: self.head_buf.as_slice().strip_prefix(b">").unwrap().strip_suffix(b"\n").unwrap(), 
                                seq: self.seq_buf.as_slice(), // Newlines are already trimmed before
                                qual: None});
        }
    }

    pub fn new(input: R, filetype: FileType) -> Self{
        FastXReader{filetype,
                    input,
                    seq_buf: Vec::<u8>::new(),
                    head_buf: Vec::<u8>::new(),
                    qual_buf: Vec::<u8>::new(),
                    plus_buf: Vec::<u8>::new(),
                    fasta_temp_buf: Vec::<u8>::new(),}
    }

    pub fn into_db(mut self) -> crate::seq_db::SeqDB{
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

        while let Some(rec) = self.read_next(){
            headbuf.extend_from_slice(rec.head);
            seqbuf.extend_from_slice(rec.seq);
            if let Some(qual) = rec.qual{
                qualbuf.as_mut().expect("Error: found a fastq record in a fasta stream.").extend_from_slice(qual);
                qual_starts.as_mut().unwrap().push(qualbuf.as_ref().unwrap().len());
            }
            head_starts.push(headbuf.len());
            seq_starts.push(seqbuf.len());
        }
        crate::seq_db::SeqDB{headbuf, seqbuf, qualbuf, head_starts, seq_starts, qual_starts}
    }



}


// Trait for a stream returning SeqRecord objects, used in DynamicFastXReader to abstract over
// The input stream type.
trait SeqRecordProducer {
    fn read_next(&mut self) -> Option<RefRecord>;

    // Since we want to call this for trait objects where we don't know the size of the struct,
    // We need to take self in a Box.
    fn into_db_boxed(self: Box<Self>) -> crate::seq_db::SeqDB;

    fn filetype(&self)-> FileType; 
}

pub struct DynamicFastXReader {
    stream: Box<dyn SeqRecordProducer>,
}

// A class that contains a dynamic trait object for different
// types of input streams.
impl DynamicFastXReader {
 
    // Need to constrain + 'static because boxed things always need to have a static
    // lifetime.
    pub fn new_from_input_stream<R: std::io::BufRead + 'static>(r: R, filetype: FileType) -> Self{
        let reader = FastXReader::<R>::new(r, filetype);
        DynamicFastXReader {stream: Box::new(reader)}
    }

    // New from file
    pub fn new_from_file(filename: &String) -> Self {
        let input = File::open(filename).unwrap();
        let (fileformat, gzipped) = figure_out_file_format(filename.as_str());
        if gzipped{

            // The GzDecoder structs have internal buffering, so we can feed in an unbuffered File stream.
            let gzdecoder = MultiGzDecoder::<File>::new(input);

            // We wrap this in BufReader because the FastX parser requires buffered reading
            let gzdecoder = BufReader::<MultiGzDecoder::<File>>::new(gzdecoder);

            Self::new_from_input_stream(gzdecoder, fileformat)

        } else{
            Self::new_from_input_stream(BufReader::<File>::new(input), fileformat)
        }
    }

    // New from stdin
    pub fn new_from_stdin(filetype: FileType, gzipped: bool) -> Self {
        let input = std::io::stdin();
        if gzipped {

            // The GzDecoder structs have internal buffering, so we can feed in an unbuffered stdin stream.
            let gzdecoder = MultiGzDecoder::<std::io::Stdin>::new(input);

            // We wrap this in BufReader because the FastX parser requires buffered reading
            let gzdecoder = BufReader::<MultiGzDecoder::<std::io::Stdin>>::new(gzdecoder);

            Self::new_from_input_stream(gzdecoder, filetype)
        } else {
            Self::new_from_input_stream(BufReader::<std::io::Stdin>::new(input), filetype)
        }
    }

    // Returns None if no more records
    pub fn read_next(&mut self) -> Option<RefRecord>{
        return self.stream.read_next()
    }

    pub fn filetype(&self)-> FileType{
        self.stream.filetype()
    }

    pub fn into_db(self) -> crate::seq_db::SeqDB{
        self.stream.into_db_boxed()
    }

}

// Implement common SeqStream trait for all
// FastXReaders over the generic parameter R.
impl<R: BufRead> SeqRecordProducer for FastXReader<R>{

    fn read_next(&mut self) -> Option<RefRecord>{
        self.read_next()
    }

    fn filetype(&self)-> FileType{
        self.filetype
    }

    fn into_db_boxed(self: Box<Self>) -> crate::seq_db::SeqDB{
        self.into_db()
    }

}

