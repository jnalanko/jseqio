
use std::{str, path::Path};

pub mod reader;
pub mod writer;
pub mod record;
pub mod seq_db;

#[derive(Copy, Clone, Debug)]
pub enum FileType{
    FASTA,
    FASTQ,
}

// Returns (file type, is_gzipped)
pub fn figure_out_file_format<P: AsRef<Path>>(filepath: P) -> (FileType, bool){
    let filename = filepath.as_ref().as_os_str().to_str().unwrap();

    let is_gzipped = filename.ends_with(".gz");
    let filename = if is_gzipped{
        &filename[0 .. filename.len()-3] // Drop the .gz suffix
    } else {filename};
    let fasta_extensions = vec![".fasta", ".fna", ".ffn", ".faa", ".frn", ".fa"];
    let fastq_extensions = vec![".fastq", ".fq"];
    if fasta_extensions.iter().any(|&suffix| filename.ends_with(suffix)){
        (FileType::FASTA, is_gzipped)
    } else if fastq_extensions.iter().any(|&suffix| filename.ends_with(suffix)){
        (FileType::FASTQ, is_gzipped)
    } else{
        panic!("Unkown file extension: {}", filename);
    }
}
