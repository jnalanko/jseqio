use std::fmt;

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

impl<'a> RefRecord<'a>{
    pub fn to_owned(&self) -> OwnedRecord{
        OwnedRecord { 
            head: self.head.to_vec(), 
            seq: self.seq.to_vec(), 
            qual: match self.qual {
                Some(q) => Some(q.to_vec()), 
                None => None
            }
        }
    }
}

impl OwnedRecord{
    pub fn as_ref_record(&self) -> RefRecord{
        RefRecord { 
            head: self.head.as_slice(), 
            seq: self.seq.as_slice(), 
            qual: match &self.qual {
                Some(q) => Some(q.as_slice()), 
                None => None
            }
        }
    }
}

impl<'a> fmt::Display for RefRecord<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f,
               "SeqRecord{{ \n  Head: {}\n  Seq:  {}\n  Qual: {}\n}}", 
               std::str::from_utf8(self.head).unwrap(),
               std::str::from_utf8(self.seq).unwrap(),
               match self.qual{
                   Some(q) => std::str::from_utf8(q).unwrap(),
                   None => "", // No quality values
               }
               
        )
    }
}