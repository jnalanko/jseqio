use std::borrow::Borrow;
use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, BufRead, BufReader};
use std::marker::PhantomData;

use flate2::read::MultiGzDecoder;

#[derive(Debug, Hash, PartialEq, Eq)]
pub struct OwnedRecord{
    pub head: Vec<u8>,    
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>, // If FASTA, this is None
}

#[derive(Debug, Hash, PartialEq, Eq, Clone, Copy)]
pub struct RefRecord<'a>{
    pub head: &'a [u8],    
    pub seq: &'a [u8],
    pub qual: Option<&'a [u8]>, // If FASTA, this is None
}

// Questions: BufRead or Read

// Functionality:
// In-memory: all().iter()
// Streaming: next()
// Get name which is the first token of the header
// Result types
// Figure out file extension string including the .gz part

struct SeqDB {
    headbuf: Vec<u8>,
    seqbuf: Vec<u8>,
    qualbuf: Vec<u8>, // Todo: option
    head_starts: Vec<usize>, // Contains end sentinel at the end
    seq_starts: Vec<usize>, // Contains end sentinel at the end
    qual_starts: Vec<usize>, // Contains end sentinel at the end
}

impl SeqDB{
    fn iter(&self) -> SeqDBIterator {
        SeqDBIterator{seq_db: self, pos: 0}
    }
}

struct SeqDBIterator<'a>{
    seq_db: &'a SeqDB,
    pos: usize,
}

impl<'a> Iterator for SeqDBIterator<'a> {
    type Item = RefRecord<'a>;

    fn next(&mut self) -> Option<RefRecord<'a>> {
        match self.pos{
            i if i < self.seq_db.head_starts.len() - 1 => {
                let head_start = self.seq_db.head_starts[i];
                let head_end = self.seq_db.head_starts[i+1];
                let seq_start = self.seq_db.seq_starts[i];
                let seq_end = self.seq_db.seq_starts[i+1];
                let qual_start = self.seq_db.qual_starts[i];
                let qual_end = self.seq_db.qual_starts[i+1];
                self.pos += 1; // Advance pointer to next element for the next round
                Some(RefRecord{
                    head: &self.seq_db.headbuf[head_start..head_end],
                    seq: &self.seq_db.seqbuf[seq_start..seq_end],
                    qual: Some(&self.seq_db.qualbuf[qual_start..qual_end]),
                })
            }
            _ => None, // End of iteration
        }
    }
}


#[derive(Copy, Clone)]
pub enum FileType{
    FASTA,
    FASTQ,
}

// Returns (file type, is_gzipped)
pub fn figure_out_file_format(filename: &str) -> (FileType, bool){
    let is_gzipped = filename.ends_with(".gz");
    let filename = if is_gzipped{
        &filename[0 .. filename.len()-3] // Drop the .gz suffix
    } else {&filename};
    let fasta_extensions = vec![".fasta", ".fna", ".ffn", ".faa", ".frn", ".fa"];
    let fastq_extensions = vec![".fastq", ".fq"];
    if fasta_extensions.iter().any(|&suffix| filename.ends_with(suffix)){
        return (FileType::FASTA, is_gzipped);
    } else if fastq_extensions.iter().any(|&suffix| filename.ends_with(suffix)){
        return (FileType::FASTQ, is_gzipped);
    } else{
        panic!("Unkown file extension: {}", filename);
    }
}

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
    pub fn next(&mut self) -> Option<RefRecord>{
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
            if self.fasta_temp_buf.len() == 0 {
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
                            if self.seq_buf.len() == 0{
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
        FastXReader{filetype: filetype,
                    input: input,
                    seq_buf: Vec::<u8>::new(),
                    head_buf: Vec::<u8>::new(),
                    qual_buf: Vec::<u8>::new(),
                    plus_buf: Vec::<u8>::new(),
                    fasta_temp_buf: Vec::<u8>::new(),}
    }

}


// Trait for a stream returning SeqRecord objects.
pub trait SeqRecordProducer {
    fn next(&mut self) -> Option<RefRecord>;
    fn filetype(&self )-> FileType; 
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
        let input = File::open(&filename).unwrap();
        let (fileformat, gzipped) = figure_out_file_format(&filename.as_str());
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
        return self.stream.next()
    }

    pub fn filetype(&self)-> FileType{
        self.stream.filetype()
    } 

}

// Implement common SeqStream trait for all
// FastXReaders over the generic parameter R.
impl<R: BufRead> SeqRecordProducer for FastXReader<R>{
    fn next(&mut self) -> Option<RefRecord>{
        self.next()
    }

    fn filetype(&self)-> FileType{
        self.filetype
    }
}

fn main() {

    //let reader = std::fs::File::open("filename.txt");
    //let buf_reader = std::io::BufReader::new(reader.unwrap());
    //let buf_buf_reader = std::io::BufReader::new(buf_reader);

/*
    let mut reader = FastXReader{
        headbuf: Vec::<u8>::new(),
        seqbuf: Vec::<u8>::new(),
        qualbuf: Vec::<u8>::new(),
        input: std::io::stdin(),
        dummy_counter: 0,
    };
    let mut records: Vec<RefRecord> = Vec::new();

    let db = reader.all();
    for rec in db.iter(){
        println!("{:?}", rec);
        records.push(rec);
    }

    for rec in db.iter(){
        println!("{:?}", rec);
    }

    dbg!(records);

*/
}
