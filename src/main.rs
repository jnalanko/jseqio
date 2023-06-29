use std::borrow::Borrow;
use std::collections::HashMap;
use std::io::Read;
use std::marker::PhantomData;

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

struct FastXReader<R: Read>{
    headbuf: Vec<u8>,
    seqbuf: Vec<u8>,
    qualbuf: Vec<u8>,
    input: R,
    dummy_counter: usize,
}

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


impl<R: Read> FastXReader<R> {
    fn all(&mut self) -> SeqDB {

        let mut headbuf = Vec::<u8>::new();
        let mut seqbuf = Vec::<u8>::new();
        let mut qualbuf = Vec::<u8>::new();

        let mut head_starts = Vec::<usize>::new();
        let mut seq_starts = Vec::<usize>::new();
        let mut qual_starts = Vec::<usize>::new();

        head_starts.push(0);
        seq_starts.push(0);
        qual_starts.push(0);
        while let Some(r) = self.next() {
            headbuf.extend(r.head);
            seqbuf.extend(r.seq);
            qualbuf.extend(r.head);
            head_starts.push(headbuf.len());
            seq_starts.push(seqbuf.len());
            qual_starts.push(qualbuf.len());
        }

        SeqDB{headbuf, seqbuf, qualbuf, head_starts, seq_starts, qual_starts}
        
        
    }

    fn next(&mut self) -> Option<RefRecord> {
        if self.dummy_counter == 3 {
            return None;
        }

        self.dummy_counter += 1;

        // "Read the data from the file"
        self.headbuf = vec![1,2,3];
        self.seqbuf = vec![1,2,3];
        self.qualbuf = vec![1,2,3];

        Some(RefRecord{
            head: &self.headbuf,
            seq: &self.seqbuf,
            qual: Some(&self.qualbuf),
        })
    }
}


fn main() {

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
}
