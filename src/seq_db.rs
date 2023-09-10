use crate::record::{RefRecord};

pub struct SeqDB {
    headbuf: Vec<u8>,
    seqbuf: Vec<u8>,
    qualbuf: Option<Vec<u8>>, // Only exists for Fastq
    head_starts: Vec<usize>, // Contains end sentinel at the end
    seq_starts: Vec<usize>, // Contains end sentinel at the end
    qual_starts: Option<Vec<usize>>, // Contains end sentinel at the end. Only exists for Fastq
}

impl SeqDB{

    pub fn iter(&self) -> SeqDBIterator {
        SeqDBIterator{seq_db: self, pos: 0}
    }

    pub fn sequence_count(&self) -> usize{
        self.head_starts.len() - 1
        // ^ The -1 is because we have an end sentinel at the end of the head_starts vector
    }

    pub fn get(&self, seq_index: usize) -> RefRecord{
        if seq_index >= self.head_starts.len(){
            panic!("SeqDB: Sequence index {} not found in database containing {} sequences", seq_index, self.sequence_count());
        }

        let head = &self.headbuf[self.head_starts[seq_index]..self.head_starts[seq_index+1]];
        let seq = &self.seqbuf[self.seq_starts[seq_index]..self.seq_starts[seq_index+1]];
        let qual = match &self.qualbuf{
            Some(buf) => { // Have quality values
                let start = self.qual_starts.as_ref().unwrap()[seq_index];   
                let end = self.qual_starts.as_ref().unwrap()[seq_index+1];
                Some(&buf[start..end])
            }
            None => None, // No quality values
        };
        RefRecord{head, seq, qual}
    }

    pub fn new_without_quality_values() -> SeqDB{
        let headbuf: Vec<u8> = Vec::new();
        let seqbuf: Vec<u8> = Vec::new();
        let head_starts: Vec<usize> = vec![0];
        let seq_starts: Vec<usize> = vec![0];
        SeqDB{headbuf, seqbuf, qualbuf: None, head_starts, seq_starts, qual_starts: None}
    }

    pub fn new() -> SeqDB{
        let mut db = Self::new_without_quality_values();
        db.qualbuf = Some(Vec::new());
        db.qual_starts = Some(vec![0]);
        db
    }

    pub fn push_record<R: crate::record::Record>(&mut self, rec: R){
        self.headbuf.extend_from_slice(rec.head());
        self.seqbuf.extend_from_slice(rec.seq());
        if self.qualbuf.is_some(){
            if let Some(qual) = rec.qual(){
                self.qualbuf.as_mut().unwrap().extend_from_slice(qual);
                self.qual_starts.as_mut().unwrap().push(self.qualbuf.as_ref().unwrap().len());
            }
        }
        self.head_starts.push(self.headbuf.len());
        self.seq_starts.push(self.seqbuf.len());
    }

    pub fn shrink_to_fit(&mut self){
        self.headbuf.shrink_to_fit();
        self.seqbuf.shrink_to_fit();
        if let Some(qualbuf) = &mut self.qualbuf{
            qualbuf.shrink_to_fit();
        }
    }
}

pub struct SeqDBIterator<'a>{
    seq_db: &'a SeqDB,
    pos: usize,
}

impl<'a> Iterator for SeqDBIterator<'a> {
    type Item = RefRecord<'a>;

    fn next(&mut self) -> Option<RefRecord<'a>> {
        match self.pos{
            i if i < self.seq_db.head_starts.len() - 1 => { // Iteration is not finished yet
                self.pos += 1; // Advance pointer to next element for the next round
                Some(self.seq_db.get(i)) // Should never be out of bounds so we unwrap the error.
            }
            _ => None, // End of iteration
        }
    }
}
