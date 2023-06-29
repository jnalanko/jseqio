








/*
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
*/
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
