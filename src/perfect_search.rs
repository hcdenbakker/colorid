extern crate bit_vec;
extern crate kmer_fa;
extern crate murmurhash64;

use bit_vec::BitVec;
use murmurhash64::murmur_hash64a;
use std;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::prelude::*;

// test with search_bigsi function show this will  introduce a ~13 second time overhead as opposed
// to writing the search out
pub fn batch_search(
    files: Vec<&str>,
    bigsi_map: &std::collections::HashMap<usize, Vec<u8>>,
    colors_accession: &std::collections::HashMap<usize, String>,
    n_ref_kmers: &std::collections::HashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
    k_size: usize,
    //filter: usize,
    cov: f64,
    gene_search: bool,
) {
    for file in files {
        //only fasta formatted file!
        eprintln!("Counting k-mers, this may take a while!");
        let vec_query = kmer_fa::read_fasta(file.to_owned());
        let kmers_query = kmer_fa::kmerize_vector(vec_query, k_size);
        //let kmers_query = kmer_fa::clean_map(unfiltered, filter);
        println!("{} kmers in query", kmers_query.len());
        let mut kmer_slices = Vec::new();
        for k in kmers_query.keys() {
            for i in 0..num_hash {
                let bit_index = murmur_hash64a(k.as_bytes(), i as u64) % bloom_size as u64;
                let bi = bit_index as usize;
                if !bigsi_map.contains_key(&bi) {
                    break;
                } else {
                    kmer_slices.push(&bigsi_map[&bi]);
                }
            }
        }
        if kmer_slices.len() < (num_hash * kmers_query.len()) {
            println!("No perfect hits!");
        } else {
            //bit-wise AND
            let mut first = BitVec::from_bytes(&kmer_slices[0].to_owned());
            for i in 1..(num_hash * kmers_query.len()) {
                let j = i as usize;
                first.intersect(&BitVec::from_bytes(&kmer_slices[j]));
            }
            let mut hits = Vec::new();
            for i in 0..first.len() {
                if first[i] {
                    hits.push(&colors_accession[&i]);
                }
            }
            println!("{} hits", hits.len());
            for h in &hits {
                println!("{}\t{}\t{}\t1.00", file, h, kmers_query.len());
            }
        }
    }
}
