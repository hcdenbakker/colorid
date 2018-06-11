extern crate bincode;
extern crate bit_vec;
extern crate flate2;
extern crate kmer_fa;
extern crate murmurhash64;
extern crate simple_bloom;

use bit_vec::BitVec;
use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use murmurhash64::murmur_hash64a;
use simple_bloom::BloomFilter;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::BufReader;
use std::io::Write;
use std::io::prelude::*;
use std::time::SystemTime;
use std;

pub fn batch_search(
    files1: Vec<&str>,
    files2: Vec<&str>,
    bigsi_map: &std::collections::HashMap<usize, Vec<u8>>,
    colors_accession: &std::collections::HashMap<usize, String>,
    n_ref_kmers: &std::collections::HashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
    k_size: usize,
    filter: isize,
    cov: f64,
    gene_search: bool,
) {
    for (i,file1) in files1.iter().enumerate() {
        if file1.ends_with("gz") {
            let unfiltered = if files2.len() == 0{
            eprintln!("{}", file1);
            eprintln!("Counting k-mers, this may take a while!");
            kmer_fa::kmers_from_fq(file1.to_owned().to_string(), k_size)
            }else{
                eprintln!("Paired end: {} {}", file1, files2[i]);
                eprintln!("Counting k-mers, this may take a while!");
                kmer_fa::kmers_fq_pe(vec![&file1, &files2[i]], k_size)
            };
            let kmers_query =
                if filter < 0{
                    let cutoff = kmer_fa::auto_cutoff(unfiltered.to_owned());
                    kmer_fa::clean_map(unfiltered, cutoff)
                }else{
                    kmer_fa::clean_map(unfiltered, filter as usize)
                };
            let num_kmers = kmers_query.len() as f64;
            println!("{} k-mers in query", num_kmers);
            let bigsi_search = SystemTime::now();
            let mut report = HashMap::new();
            let mut uniq_freqs = HashMap::new();
            for k in kmers_query.keys() {
                let mut kmer_slices = Vec::new();
                for i in 0..num_hash {
                    let bit_index = murmur_hash64a(k.as_bytes(), i as u64) % bloom_size as u64;
                    let bi = bit_index as usize;
                    if !bigsi_map.contains_key(&bi) {
                        break;
                    } else {
                        kmer_slices.push(&bigsi_map[&bi]);
                    }
                }
                if kmer_slices.len() < num_hash {
                    continue;
                } else {
                    let mut first = BitVec::from_bytes(&kmer_slices[0].to_owned());
                    for i in 1..num_hash {
                        let j = i as usize;
                        first.intersect(&BitVec::from_bytes(&kmer_slices[j]));
                    }
                    let mut hits = Vec::new();
                    for i in 0..first.len() {
                        if first[i] {
                            hits.push(&colors_accession[&i]);
                        }
                    }
                    for h in &hits {
                        let count = report.entry(h.to_string()).or_insert(0);
                        *count += 1;
                    }
                    if hits.len() == 1 {
                        let key = hits[0];
                        let value = kmers_query[&k.to_string()] as f64;
                        uniq_freqs
                            .entry(key.to_string())
                            .or_insert_with(Vec::new)
                            .push(value);
                    }
                }
            }
            match bigsi_search.elapsed() {
                Ok(elapsed) => {
                    eprintln!("Search: {} sec", elapsed.as_secs());
                }
                Err(e) => {
                    // an error occurred!
                    eprintln!("Error: {:?}", e);
                }
            }
            if !gene_search {
                super::generate_report(file1, report, &uniq_freqs, &n_ref_kmers, cov);
            } else {
                super::generate_report_gene(file1, report, num_kmers as usize);
            }
        } else {
            //else we assume it is a fasta formatted file!
            eprintln!("{}", file1);
            eprintln!("Counting k-mers, this may take a while!");
            let vec_query = kmer_fa::read_fasta(file1.to_owned().to_string());
            let unfiltered = kmer_fa::kmerize_vector(vec_query, k_size);
            let kmers_query =
                if gene_search{
                    kmer_fa::clean_map(unfiltered, 0)
                } else if filter < 0{
                    eprintln!("no gene search");
                    let cutoff = kmer_fa::auto_cutoff(unfiltered.to_owned());
                    kmer_fa::clean_map(unfiltered, cutoff)
                }else{
                    kmer_fa::clean_map(unfiltered, filter as usize)
                };
            let num_kmers = kmers_query.len() as f64;
            eprintln!("{} k-mers in query", num_kmers);
            let mut report = HashMap::new();
            let mut uniq_freqs = HashMap::new();
            for k in kmers_query.keys() {
                let mut kmer_slices = Vec::new();
                for i in 0..num_hash {
                    let bit_index = murmur_hash64a(k.as_bytes(), i as u64) % bloom_size as u64;
                    let bi = bit_index as usize;
                    if !bigsi_map.contains_key(&bi) {
                        break;
                    } else {
                        kmer_slices.push(&bigsi_map[&bi]);
                    }
                }
                if kmer_slices.len() < num_hash {
                    continue;
                } else {
                    let mut first = BitVec::from_bytes(&kmer_slices[0].to_owned());
                    for i in 1..num_hash {
                        let j = i as usize;
                        first.intersect(&BitVec::from_bytes(&kmer_slices[j]));
                    }
                    let mut hits = Vec::new();
                    for i in 0..first.len() {
                        if first[i] {
                            hits.push(&colors_accession[&i]);
                        }
                    }
                    for h in &hits {
                        let count = report.entry(h.to_string()).or_insert(0);
                        *count += 1;
                    }
                    if hits.len() == 1 {
                        let key = hits[0];
                        let value = kmers_query[&k.to_string()] as f64;
                        uniq_freqs
                            .entry(key.to_string())
                            .or_insert_with(Vec::new)
                            .push(value);
                    }
                }
            }
            if !gene_search {
                super::generate_report(file1, report, &uniq_freqs, &n_ref_kmers, cov);
            } else {
                super::generate_report_gene(file1, report, num_kmers as usize);
            }
        }
    }
}
