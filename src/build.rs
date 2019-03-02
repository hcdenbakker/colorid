extern crate bit_vec;
extern crate rayon;

use bit_vec::BitVec;
use fnv;
use kmer;
use rayon::prelude::*;
use simple_bloom;
use std;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::prelude::*;

pub fn tab_to_map(filename: String) -> fnv::FnvHashMap<std::string::String, Vec<String>> {
    let mut map = fnv::FnvHashMap::default();
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
    map: &fnv::FnvHashMap<std::string::String, Vec<String>>,
    bloom_size: usize,
    num_hash: usize,
    k_size: usize,
    quality: u8,
    cutoff: isize,
) -> (
    fnv::FnvHashMap<usize, BitVec>,
    fnv::FnvHashMap<usize, String>,
    fnv::FnvHashMap<String, usize>,
) {
    //get a hash map with taxa and bit vector from bloom filter
    let map_length = map.len();
    let mut bit_map = HashMap::with_capacity(map_length);
    let mut ref_kmer = fnv::FnvHashMap::default();
    let mut accessions = Vec::with_capacity(map_length);
    let mut counter = 1;
    for (accession, v) in map {
        eprintln!("Adding {} to index ({}/{})", accession, counter, map_length);
        counter += 1;
        if v.len() == 2 {
            let unfiltered = kmer::kmers_fq_pe_qual(vec![&v[0], &v[1]], k_size, 1, quality);
            let kmers = if cutoff == -1 {
                let auto_cutoff = kmer::auto_cutoff(unfiltered.to_owned());
                kmer::clean_map(unfiltered, auto_cutoff)
            } else {
                kmer::clean_map(unfiltered, cutoff as usize)
            };
            ref_kmer.insert(accession.to_string(), kmers.len());
            let mut filter = simple_bloom::BloomFilter::new(bloom_size as usize, num_hash as usize);
            for kmer in kmers.keys() {
                filter.insert(&kmer);
            }
            bit_map.insert(accession, filter.bits);
        } else {
            if v[0].ends_with("gz") {
                let unfiltered = kmer::kmers_from_fq_qual(v[0].to_owned(), k_size, 1, quality);
                let kmers = if cutoff == -1 {
                    let auto_cutoff = kmer::auto_cutoff(unfiltered.to_owned());
                    kmer::clean_map(unfiltered, auto_cutoff)
                } else {
                    kmer::clean_map(unfiltered, cutoff as usize)
                };
                ref_kmer.insert(accession.to_string(), kmers.len());
                let mut filter =
                    simple_bloom::BloomFilter::new(bloom_size as usize, num_hash as usize);
                for kmer in kmers.keys() {
                    filter.insert(&kmer);
                }
                bit_map.insert(accession, filter.bits);
            } else {
                let vec = kmer::read_fasta(v[0].to_string());
                let kmers = if cutoff == -1 {
                    kmer::kmerize_vector(&vec, k_size, 1)
                } else {
                    let unfiltered = kmer::kmerize_vector(&vec, k_size, 1);
                    kmer::clean_map(unfiltered, cutoff as usize)
                };
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
    accessions.sort();
    //create hash table with colors for accessions
    let num_taxa = accessions.len();
    let mut accession_colors = fnv::FnvHashMap::default();
    let mut colors_accession = fnv::FnvHashMap::default();
    for (c, s) in accessions.iter().enumerate() {
        accession_colors.insert(s.to_string(), c);
        colors_accession.insert(c, s.to_string());
    }
    //create actual index, the most straight forward way, but not very efficient
    let mut bigsi_map = fnv::FnvHashMap::default();
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
            bigsi_map.insert(i, bitvec);
        }
    }
    (bigsi_map, colors_accession, ref_kmer)
}

pub fn build_multi(
    map: &fnv::FnvHashMap<std::string::String, Vec<String>>,
    bloom_size: usize,
    num_hash: usize,
    k_size: usize,
    t: usize,
    quality: u8,
    cutoff: isize,
) -> (
    fnv::FnvHashMap<usize, BitVec>,
    fnv::FnvHashMap<usize, String>,
    fnv::FnvHashMap<String, usize>,
) {
    //get a hash map with taxa and bit vector from bloom filter
    rayon::ThreadPoolBuilder::new()
        .num_threads(t)
        .build_global()
        .unwrap();
    let map_length = map.len();
    let mut ref_kmer = fnv::FnvHashMap::default();
    let mut accessions = Vec::with_capacity(map_length);
    let mut map_vec = Vec::with_capacity(map_length);
    for (accession, v) in map {
        map_vec.push((accession, v));
    }
    let mut colors_accession = fnv::FnvHashMap::default();
    let mut accession_colors = fnv::FnvHashMap::default();
    let bloom_vec: Vec<usize> = (0..bloom_size).collect();
    let d: Vec<_>;
    {
        let c: Vec<_>;
        eprintln!(
            "Inference of Bloom filters in parallel using {} threads.",
            t
        );
        c = map_vec
            .par_iter()
            .map(|l| {
                if l.1.len() == 2 {
                    let unfiltered =
                        kmer::kmers_fq_pe_qual(vec![&l.1[0], &l.1[1]], k_size, 1, quality);
                    let kmers = if cutoff == -1 {
                        let auto_cutoff = kmer::auto_cutoff(unfiltered.to_owned());
                        kmer::clean_map(unfiltered, auto_cutoff)
                    } else {
                        kmer::clean_map(unfiltered, cutoff as usize)
                    };
                    let mut filter =
                        simple_bloom::BloomFilter::new(bloom_size as usize, num_hash as usize);
                    for kmer in kmers.keys() {
                        filter.insert(&kmer);
                    }
                    (l.0, filter.bits, kmers.len())
                } else {
                    if l.1[0].ends_with("gz") {
                        let unfiltered = kmer::kmers_from_fq_qual(l.1[0].to_owned(), k_size, 1, 15);
                        let kmers = if cutoff == -1 {
                            let auto_cutoff = kmer::auto_cutoff(unfiltered.to_owned());
                            kmer::clean_map(unfiltered, auto_cutoff)
                        } else {
                            kmer::clean_map(unfiltered, cutoff as usize)
                        };
                        let mut filter =
                            simple_bloom::BloomFilter::new(bloom_size as usize, num_hash as usize);
                        for kmer in kmers.keys() {
                            filter.insert(&kmer);
                        }
                        (l.0, filter.bits, kmers.len())
                    } else {
                        let vec = kmer::read_fasta(l.1[0].to_string());
                        let kmers = if cutoff == -1 {
                            kmer::kmerize_vector(&vec, k_size, 1)
                        } else {
                            let unfiltered = kmer::kmerize_vector(&vec, k_size, 1);
                            kmer::clean_map(unfiltered, cutoff as usize)
                        };
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
        for t in &c {
            ref_kmer.insert(t.0.to_string(), t.2);
        }
        for accession in map.keys() {
            accessions.push(accession);
        }
        accessions.sort();
        //create hash table with colors for accessions
        let num_taxa = accessions.len();
        for (c, s) in accessions.iter().enumerate() {
            accession_colors.insert(s.to_string(), c);
            colors_accession.insert(c, s.to_string());
        }
        eprintln!("Creation of index, this may take a while!");
        //create actual index, the most straight forward way, but not very efficient
        d = bloom_vec
            .par_iter()
            .map(|i| {
                let mut bitvec = BitVec::from_elem(num_taxa, false);
                for a in &c {
                    if a.1[*i] {
                        bitvec.set(accession_colors[&a.0.to_string()], true);
                    }
                }
                (i, bitvec)
            })
            .collect();
        drop(c);
    }
    let mut bigsi_map = fnv::FnvHashMap::default();
    for t in d {
        if t.1.none() {
            continue;
        } else {
            bigsi_map.insert(t.0.to_owned(), t.1);
        }
    }
    (bigsi_map, colors_accession, ref_kmer)
}

pub fn build_multi_mini(
    map: &fnv::FnvHashMap<std::string::String, Vec<String>>,
    bloom_size: usize,
    num_hash: usize,
    k_size: usize,
    m: usize,
    t: usize,
    quality: u8,
    cutoff: isize,
) -> (
    fnv::FnvHashMap<usize, BitVec>,
    fnv::FnvHashMap<usize, String>,
    fnv::FnvHashMap<String, usize>,
) {
    //get a hash map with taxa and bit vector from bloom filter
    rayon::ThreadPoolBuilder::new()
        .num_threads(t)
        .build_global()
        .unwrap();
    let map_length = map.len();
    let mut ref_kmer = fnv::FnvHashMap::default();
    let mut accessions = Vec::with_capacity(map_length);
    let mut map_vec = Vec::with_capacity(map_length);
    for (accession, v) in map {
        map_vec.push((accession, v));
    }
    let mut colors_accession = fnv::FnvHashMap::default();
    let mut accession_colors = fnv::FnvHashMap::default();
    let bloom_vec: Vec<usize> = (0..bloom_size).collect();
    let d: Vec<_>;
    {
        let c: Vec<_>;
        eprintln!(
            "Inference of Bloom filters in parallel using {} threads.",
            t
        );
        c = map_vec
            .par_iter()
            .map(|l| {
                if l.1.len() == 2 {
                    let unfiltered = kmer::kmers_fq_pe_minimizer_qual(
                        vec![&l.1[0], &l.1[1]],
                        k_size,
                        m,
                        1,
                        quality,
                    );
                    let kmers = if cutoff == -1 {
                        let auto_cutoff = kmer::auto_cutoff(unfiltered.to_owned());
                        kmer::clean_map(unfiltered, auto_cutoff)
                    } else {
                        kmer::clean_map(unfiltered, cutoff as usize)
                    };
                    let mut filter =
                        simple_bloom::BloomFilter::new(bloom_size as usize, num_hash as usize);
                    for kmer in kmers.keys() {
                        filter.insert(&kmer);
                    }
                    (l.0, filter.bits, kmers.len())
                } else {
                    if l.1[0].ends_with("gz") {
                        let unfiltered = kmer::kmers_from_fq_minimizer_qual(
                            l.1[0].to_owned(),
                            k_size,
                            m,
                            1,
                            quality,
                        );
                        let kmers = if cutoff == -1 {
                            let auto_cutoff = kmer::auto_cutoff(unfiltered.to_owned());
                            kmer::clean_map(unfiltered, auto_cutoff)
                        } else {
                            kmer::clean_map(unfiltered, cutoff as usize)
                        };
                        let mut filter =
                            simple_bloom::BloomFilter::new(bloom_size as usize, num_hash as usize);
                        for kmer in kmers.keys() {
                            filter.insert(&kmer);
                        }
                        (l.0, filter.bits, kmers.len())
                    } else {
                        let vec = kmer::read_fasta(l.1[0].to_string());
                        let kmers = if cutoff == -1 {
                            kmer::minimerize_vector_skip_n(&vec, k_size, m, 1)
                        } else {
                            let unfiltered = kmer::minimerize_vector_skip_n(&vec, k_size, m, 1);
                            kmer::clean_map(unfiltered, cutoff as usize)
                        };
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
        for t in &c {
            ref_kmer.insert(t.0.to_string(), t.2);
        }
        for accession in map.keys() {
            accessions.push(accession);
        }
        accessions.sort();
        //create hash table with colors for accessions
        let num_taxa = accessions.len();
        for (c, s) in accessions.iter().enumerate() {
            accession_colors.insert(s.to_string(), c);
            colors_accession.insert(c, s.to_string());
        }
        eprintln!("Creation of index, this may take a while!");
        //create actual index, the most straight forward way, but not very efficient
        //this can be done in parallel!
        d = bloom_vec
            .par_iter()
            .map(|i| {
                let mut bitvec = BitVec::from_elem(num_taxa, false);
                for a in &c {
                    if a.1[*i] {
                        bitvec.set(accession_colors[&a.0.to_string()], true);
                    }
                }
                (i, bitvec)
            })
            .collect();
    }
    let mut bigsi_map = fnv::FnvHashMap::default();
    for t in d {
        if t.1.none() {
            continue;
        } else {
            bigsi_map.insert(t.0.to_owned(), t.1);
        }
    }
    (bigsi_map, colors_accession, ref_kmer)
}

pub fn build_single_mini(
    map: &fnv::FnvHashMap<std::string::String, Vec<String>>,
    bloom_size: usize,
    num_hash: usize,
    k_size: usize,
    m: usize,
    quality: u8,
    cutoff: isize,
) -> (
    fnv::FnvHashMap<usize, BitVec>,
    fnv::FnvHashMap<usize, String>,
    fnv::FnvHashMap<String, usize>,
) {
    //get a hash map with taxa and bit vector from bloom filter
    let map_length = map.len();
    let mut bit_map = HashMap::with_capacity(map_length);
    let mut ref_kmer = fnv::FnvHashMap::default();
    let mut accessions = Vec::with_capacity(map_length);
    let mut counter = 1;
    for (accession, v) in map {
        eprintln!("Adding {} to index ({}/{})", accession, counter, map_length);
        counter += 1;
        if v.len() == 2 {
            let unfiltered = kmer::kmers_fq_pe_qual(vec![&v[0], &v[1]], k_size, 1, quality);
            let kmers = if cutoff == -1 {
                let auto_cutoff = kmer::auto_cutoff(unfiltered.to_owned());
                kmer::clean_map(unfiltered, auto_cutoff)
            } else {
                kmer::clean_map(unfiltered, cutoff as usize)
            };
            let mut filter = simple_bloom::BloomFilter::new(bloom_size as usize, num_hash as usize);
            for kmer in kmers.keys() {
                filter.insert(&kmer::find_minimizer(&kmer, m));
            }
            bit_map.insert(accession, filter.bits);
        } else {
            if v[0].ends_with("gz") {
                let unfiltered = kmer::kmers_from_fq_qual(v[0].to_owned(), k_size, 1, quality);
                let kmers = if cutoff == -1 {
                    let auto_cutoff = kmer::auto_cutoff(unfiltered.to_owned());
                    kmer::clean_map(unfiltered, auto_cutoff)
                } else {
                    kmer::clean_map(unfiltered, cutoff as usize)
                };
                let mut filter =
                    simple_bloom::BloomFilter::new(bloom_size as usize, num_hash as usize);
                for kmer in kmers.keys() {
                    filter.insert(&kmer::find_minimizer(&kmer, m));
                }
                bit_map.insert(accession, filter.bits);
            } else {
                let vec = kmer::read_fasta(v[0].to_string());
                let kmers = if cutoff == -1 {
                    kmer::kmerize_vector(&vec, k_size, 1)
                } else {
                    let unfiltered = kmer::kmerize_vector(&vec, k_size, 1);
                    kmer::clean_map(unfiltered, cutoff as usize)
                };
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
    accessions.sort();
    //create hash table with colors for accessions
    let num_taxa = accessions.len();
    let mut accession_colors = fnv::FnvHashMap::default();
    let mut colors_accession = fnv::FnvHashMap::default();
    for (c, s) in accessions.iter().enumerate() {
        accession_colors.insert(s.to_string(), c);
        colors_accession.insert(c, s.to_string());
    }
    //create actual index, the most straight forward way, but not very efficient
    let mut bigsi_map = fnv::FnvHashMap::default();
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
            bigsi_map.insert(i, bitvec);
        }
    }
    (bigsi_map, colors_accession, ref_kmer)
}
