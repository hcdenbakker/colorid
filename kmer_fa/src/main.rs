extern crate kmer_fa;

use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    let filename = &args[1];
    if filename.ends_with("gz") {
        let map = kmer_fa::kmers_from_fq(filename.to_string(), 27);
        println!("{}", map.len());
        let clean = kmer_fa::clean_map(map, 3);
        println!("{}", clean.len());
    } else {
        let vec = kmer_fa::read_fasta(filename.to_string());
        let map = kmer_fa::kmerize_vector(vec, 27);
        println!("{}", map.len())
    }
}
