extern crate rayon;

use bit_vec::BitVec;
use kmer;
use rayon::prelude::*;
use simple_bloom;
use std;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::prelude::*;

pub fn tab_to_map(filename: String) -> std::collections::HashMap<std::string::String, Vec<String>> {
    let mut map = HashMap::new();
    let f = File::open(filename).expect("reference file not found");
    for line in io::BufReader::new(f).lines() {
        let l = line.unwrap();
        let v: Vec<&str> = l.split('\t').collect();
        if v.len() == 2 {
            map.insert(String::from(v[0]), vec![String::from(v[1])]);
        } else {
            map.insert(
                String::from(v[0]),
                vec![String::from(v[1]), String::from(v[2])],
            );
        }
    }
    map
}

pub fn build_single(
    map: &std::collections::HashMap<std::string::String, Vec<String>>,
    bloom_size: usize,
    num_hash: usize,
    k_size: usize,
) -> (
    //std::collections::HashMap<usize, bit_vec::BitVec>,
    std::collections::HashMap<usize, Vec<u8>>,
    std::collections::HashMap<usize, String>,
    std::collections::HashMap<String, usize>,
) {
    //get a hash map with taxa and bit vector from bloom filter
    let map_length = map.len();
    let mut bit_map = HashMap::with_capacity(map_length);
    let mut ref_kmer = HashMap::with_capacity(map_length);
    let mut accessions = Vec::with_capacity(map_length);
    let mut counter = 1;
    for (accession, v) in map {
        eprintln!("Adding {} to index ({}/{})", accession, counter, map_length);
        counter += 1;
        if v.len() == 2 {
            let unfiltered = kmer::kmers_fq_pe(vec![&v[0], &v[1]], k_size);
            let cutoff = kmer::auto_cutoff(unfiltered.to_owned());
            let kmers = kmer::clean_map(unfiltered, cutoff);
            let mut filter = simple_bloom::BloomFilter::new(bloom_size as usize, num_hash as usize);
            for kmer in kmers.keys() {
                filter.insert(&kmer);
            }
            bit_map.insert(accession, filter.bits);
        } else {
            if v[0].ends_with("gz") {
                let unfiltered = kmer::kmers_from_fq(v[0].to_owned(), k_size);
                let kmers = kmer::clean_map(unfiltered, 1);
                let mut filter =
                    simple_bloom::BloomFilter::new(bloom_size as usize, num_hash as usize);
                for kmer in kmers.keys() {
                    filter.insert(&kmer);
                }
                bit_map.insert(accession, filter.bits);
            } else {
                let vec = kmer::read_fasta(v[0].to_string());
                let kmers = kmer::kmerize_vector(vec, k_size);
                ref_kmer.insert(accession.to_string(), kmers.len());
                let mut filter =
                    simple_bloom::BloomFilter::new(bloom_size as usize, num_hash as usize);
                for kmer in kmers.keys() {
                    filter.insert(&kmer);
                }
                bit_map.insert(accession, filter.bits);
            }
        }
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
    //create actual index, the most straight forward way, but not very efficient
    let mut bigsi_map = HashMap::new();
    for i in 0..bloom_size {
        let mut bitvec = BitVec::from_elem(num_taxa, false);
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

pub fn build_multi(
    map: &std::collections::HashMap<std::string::String, Vec<String>>,
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
            if l.1.len() == 2 {
                let unfiltered = kmer::kmers_fq_pe(vec![&l.1[0], &l.1[1]], k_size);
                let cutoff = kmer::auto_cutoff(unfiltered.to_owned());
                let kmers = kmer::clean_map(unfiltered, cutoff);
                let mut filter =
                    simple_bloom::BloomFilter::new(bloom_size as usize, num_hash as usize);
                for kmer in kmers.keys() {
                    filter.insert(&kmer);
                }
                (l.0, filter.bits, kmers.len())
            } else {
                if l.1[0].ends_with("gz") {
                    let unfiltered = kmer::kmers_from_fq(l.1[0].to_owned(), k_size);
                    let cutoff = kmer::auto_cutoff(unfiltered.to_owned());
                    let kmers = kmer::clean_map(unfiltered, cutoff);
                    let mut filter =
                        simple_bloom::BloomFilter::new(bloom_size as usize, num_hash as usize);
                    for kmer in kmers.keys() {
                        filter.insert(&kmer);
                    }
                    (l.0, filter.bits, kmers.len())
                } else {
                    let vec = kmer::read_fasta(l.1[0].to_string());
                    let kmers = kmer::kmerize_vector(vec, k_size);
                    let mut filter =
                        simple_bloom::BloomFilter::new(bloom_size as usize, num_hash as usize);
                    for kmer in kmers.keys() {
                        filter.insert(&kmer);
                    }
                    (l.0, filter.bits, kmers.len())
                }
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
        let mut bitvec = BitVec::from_elem(num_taxa, false);
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

pub fn build_multi_mini(
    map: &std::collections::HashMap<std::string::String, Vec<String>>,
    bloom_size: usize,
    num_hash: usize,
    k_size: usize,
    m: usize,
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
            if l.1.len() == 2 {
                let unfiltered = kmer::kmers_fq_pe_minimizer(vec![&l.1[0], &l.1[1]], k_size, m);
                let cutoff = kmer::auto_cutoff(unfiltered.to_owned());
                let kmers = kmer::clean_map(unfiltered, cutoff);
                let mut filter =
                    simple_bloom::BloomFilter::new(bloom_size as usize, num_hash as usize);
                for kmer in kmers.keys() {
                    filter.insert(&kmer);
                }
                (l.0, filter.bits, kmers.len())
            } else {
                if l.1[0].ends_with("gz") {
                    let unfiltered = kmer::kmers_from_fq_minimizer(l.1[0].to_owned(), k_size, m);
                    let cutoff = kmer::auto_cutoff(unfiltered.to_owned());
                    let kmers = kmer::clean_map(unfiltered, cutoff);
                    let mut filter =
                        simple_bloom::BloomFilter::new(bloom_size as usize, num_hash as usize);
                    for kmer in kmers.keys() {
                        filter.insert(&kmer);
                    }
                    (l.0, filter.bits, kmers.len())
                } else {
                    let vec = kmer::read_fasta(l.1[0].to_string());
                    let kmers = kmer::minimerize_vector(vec, k_size, m);
                    let mut filter =
                        simple_bloom::BloomFilter::new(bloom_size as usize, num_hash as usize);
                    for kmer in kmers.keys() {
                        filter.insert(&kmer);
                    }
                    (l.0, filter.bits, kmers.len())
                }
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
        let mut bitvec = BitVec::from_elem(num_taxa, false);
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

pub fn build_single_mini(
    map: &std::collections::HashMap<std::string::String, Vec<String>>,
    bloom_size: usize,
    num_hash: usize,
    k_size: usize,
    m: usize,
) -> (
    //std::collections::HashMap<usize, bit_vec::BitVec>,
    std::collections::HashMap<usize, Vec<u8>>,
    std::collections::HashMap<usize, String>,
    std::collections::HashMap<String, usize>,
) {
    //get a hash map with taxa and bit vector from bloom filter
    let map_length = map.len();
    let mut bit_map = HashMap::with_capacity(map_length);
    let mut ref_kmer = HashMap::with_capacity(map_length);
    let mut accessions = Vec::with_capacity(map_length);
    let mut counter = 1;
    for (accession, v) in map {
        eprintln!("Adding {} to index ({}/{})", accession, counter, map_length);
        counter += 1;
        if v.len() == 2 {
            let unfiltered = kmer::kmers_fq_pe(vec![&v[0], &v[1]], k_size);
            let cutoff = kmer::auto_cutoff(unfiltered.to_owned());
            let kmers = kmer::clean_map(unfiltered, cutoff);
            let mut filter = simple_bloom::BloomFilter::new(bloom_size as usize, num_hash as usize);
            for kmer in kmers.keys() {
                filter.insert(&kmer::find_minimizer(&kmer, m));
            }
            bit_map.insert(accession, filter.bits);
        } else {
            if v[0].ends_with("gz") {
                let unfiltered = kmer::kmers_from_fq(v[0].to_owned(), k_size);
                let kmers = kmer::clean_map(unfiltered, 1);
                let mut filter =
                    simple_bloom::BloomFilter::new(bloom_size as usize, num_hash as usize);
                for kmer in kmers.keys() {
                    filter.insert(&kmer::find_minimizer(&kmer, m));
                }
                bit_map.insert(accession, filter.bits);
            } else {
                let vec = kmer::read_fasta(v[0].to_string());
                let kmers = kmer::kmerize_vector(vec, k_size);
                ref_kmer.insert(accession.to_string(), kmers.len());
                let mut filter =
                    simple_bloom::BloomFilter::new(bloom_size as usize, num_hash as usize);
                for kmer in kmers.keys() {
                    filter.insert(&kmer::find_minimizer(&kmer, m));
                }
                bit_map.insert(accession, filter.bits);
            }
        }
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
    //create actual index, the most straight forward way, but not very efficient
    let mut bigsi_map = HashMap::new();
    for i in 0..bloom_size {
        let mut bitvec = BitVec::from_elem(num_taxa, false);
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
