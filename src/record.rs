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
