use crate::record::{RefRecord};

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

    pub fn get(&self, seq_index: usize) -> Result<RefRecord, Box<dyn std::error::Error>>{
        if seq_index >= self.head_starts.len(){
            let msg = format!("SeqDB: Sequence index {} not found in database containing {} sequences", seq_index, self.head_starts.len()-1);
            // ^ The -1 is because we have an end sentinel at the end of the head_starts vector
            return Err(msg.into());
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
        Ok(RefRecord{head, seq, qual})
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
                Some(self.seq_db.get(i).unwrap()) // Should never be out of bounds so we unwrap the error.
            }
            _ => None, // End of iteration
        }
    }
}
