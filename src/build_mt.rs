extern crate bit_vec;
extern crate flate2;
extern crate kmer_fa;
extern crate murmurhash64;
extern crate rayon;
extern crate simple_bloom;

use bit_vec::BitVec;
use flate2::read::MultiGzDecoder;
use murmurhash64::murmur_hash64a;
use rayon::prelude::*;
use simple_bloom::BloomFilter;
use std;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::prelude::*;

pub fn build_bigsi(
    map: &std::collections::HashMap<std::string::String, String>,
    bloom_size: usize,
    num_hash: usize,
    k_size: usize,
    t: usize,
) -> (
    //std::collections::HashMap<usize, bit_vec::BitVec>,
    std::collections::HashMap<usize, Vec<u8>>,
    std::collections::HashMap<usize, String>,
    std::collections::HashMap<String, usize>,
) {
    //get a hash map with taxa and bit vector from bloom filter
    rayon::ThreadPoolBuilder::new()
        .num_threads(t)
        .build_global()
        .unwrap();
    let map_length = map.len();
    let mut bit_map = HashMap::with_capacity(map_length);
    let mut ref_kmer = HashMap::with_capacity(map_length);
    let mut accessions = Vec::with_capacity(map_length);
    let mut map_vec = Vec::with_capacity(map_length);
    for (accession, v) in map {
        map_vec.push((accession, v));
    }
    let c: Vec<_>;
    eprintln!(
        "Inference of Bloom filters in parallel using {} threads.",
        t
    );
    c = map_vec
        .par_iter()
        .map(|l| {
            if l.1.ends_with("gz") {
                let unfiltered = kmer_fa::kmers_from_fq(l.1.to_owned(), k_size);
                let kmers = kmer_fa::clean_map(unfiltered, 1);
                let mut filter = BloomFilter::new(bloom_size as usize, num_hash as usize);
                for kmer in kmers.keys() {
                    filter.insert(&kmer);
                }
                (l.0, filter.bits, kmers.len())
            } else {
                let vec = kmer_fa::read_fasta(l.1.to_string());
                let kmers = kmer_fa::kmerize_vector(vec, k_size);
                let mut filter = BloomFilter::new(bloom_size as usize, num_hash as usize);
                for kmer in kmers.keys() {
                    filter.insert(&kmer);
                }
                (l.0, filter.bits, kmers.len())
            }
        })
        .collect();
    for t in c {
        bit_map.insert(t.0, t.1);
        ref_kmer.insert(t.0.to_string(), t.2);
    }
    for accession in map.keys() {
        accessions.push(accession);
    }
    //create hash table with colors for accessions
    let num_taxa = accessions.len();
    let mut accession_colors = HashMap::with_capacity(accessions.len());
    let mut colors_accession = HashMap::with_capacity(accessions.len());
    for (c, s) in accessions.iter().enumerate() {
        accession_colors.insert(s.to_string(), c);
        colors_accession.insert(c, s.to_string());
    }
    eprintln!("Creation of index, this may take a while!");
    //create actual index, the most straight forward way, but not very efficient
    let mut bigsi_map = HashMap::new();
    for i in 0..bloom_size {
        let mut bitvec = bit_vec::BitVec::from_elem(num_taxa, false);
        for (t, s) in &bit_map {
            if s[i] {
                bitvec.set(accession_colors[&t.to_string()], true);
            }
        }
        if bitvec.none() {
            continue;
        } else {
            bigsi_map.insert(i, bitvec.to_bytes());
        }
    }
    (bigsi_map, colors_accession, ref_kmer)
}
