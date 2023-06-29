use crate::record::{RefRecord, Record};

pub struct SeqDB {
    pub headbuf: Vec<u8>,
    pub seqbuf: Vec<u8>,
    pub qualbuf: Option<Vec<u8>>, // Only exists for Fastq
    pub head_starts: Vec<usize>, // Contains end sentinel at the end
    pub seq_starts: Vec<usize>, // Contains end sentinel at the end
    pub qual_starts: Option<Vec<usize>>, // Contains end sentinel at the end. Only exists for Fastq
}

impl SeqDB{
    pub fn iter(&self) -> SeqDBIterator {
        SeqDBIterator{seq_db: self, pos: 0}
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
                let head = &self.seq_db.headbuf[self.seq_db.head_starts[i]..self.seq_db.head_starts[i+1]];
                let seq = &self.seq_db.seqbuf[self.seq_db.seq_starts[i]..self.seq_db.seq_starts[i+1]];
                let qual = match &self.seq_db.qualbuf{
                    Some(buf) => { // Have quality values
                        let start = self.seq_db.qual_starts.as_ref().unwrap()[i];   
                        let end = self.seq_db.qual_starts.as_ref().unwrap()[i+1];
                        Some(&buf[start..end])
                    }
                    None => None, // No quality values
                };

                self.pos += 1; // Advance pointer to next element for the next round
                Some(RefRecord{head, seq, qual})
            }
            _ => None, // End of iteration
        }
    }
}
