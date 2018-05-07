extern crate flate2;
extern crate murmurhash64;
extern crate kmer_fa;
extern crate bit_vec;

use std;
use std::io;
use std::io::prelude::*;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::collections::HashMap;
use murmurhash64::murmur_hash64a;
use bit_vec::BitVec;

// test with search_bigsi function show this will  introduce a ~13 second time overhead as opposed
// to writing the search out
pub fn batch_search(
    files: Vec<&str>,
    bigsi_map: std::collections::HashMap<usize, Vec<u8>>,
    colors_accession: std::collections::HashMap<usize, String>,
    n_ref_kmers: std::collections::HashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
    k_size: usize,
    filter: i32,
    cov: f64,
    gene_search: bool,
) {
    for file in files {
        if file.ends_with("gz") {
            println!("{}", file);
            eprintln!("Counting k-mers, this may take a while!");
            let unfiltered = kmer_fa::kmers_from_fq(file.to_owned(), k_size);
            let kmers_query = kmer_fa::clean_map(unfiltered, filter);
            let num_kmers = kmers_query.len() as f64;
            println!("{} k-mers in query", num_kmers);
            let bigsi_search = super::SystemTime::now();
            let mut pre_report: std::collections::HashMap<usize, usize> = HashMap::new();
            //let mut pre_report = HashMap::new();
            let mut no_hits = 0;
            let mut uniq_freqs = HashMap::new();
            for (k, _) in &kmers_query {
                let mut kmer_slices = Vec::new();
                for i in 0..num_hash {
                    let bit_index = murmur_hash64a(k.as_bytes(), i as u64) % bloom_size as u64;
                    let bi = bit_index as usize;
                    if bigsi_map.contains_key(&bi) == false {
                        no_hits += 1;
                        break;
                    } else {
                        kmer_slices.push(bigsi_map.get(&bi).unwrap());
                    }
                }
                let mut first = BitVec::from_bytes(&kmer_slices[0].to_owned());
                for i in 1..num_hash {
                    let j = i as usize;
                    first.intersect(&BitVec::from_bytes(&kmer_slices[j]));
                }
                let mut color = 0;
                let mut hits = Vec::new();
                for i in first{
                    if i == true {
                        hits.push(color);
                    }
                    color +=1;
                }
                for h in &hits {
                    *pre_report.entry(*h).or_insert(0) += 1; 
                    //let count = pre_report.entry(*h).get().unwrap_or_else(|v| v.insert(0));
                    //let count = pre_report.entry(h.or_insert(0));
                    //*count += 1;
                }
                if hits.len() == 1 {
                    let key = colors_accession.get(&hits[0]).unwrap().to_string();
                    let value = *kmers_query.get(&k.to_string()).unwrap() as f64;
                    uniq_freqs
                        .entry(key.to_string())
                        .or_insert(Vec::new())
                        .push(value);
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
            let mut report = HashMap::new();
            report.insert("no_hits".to_string(), no_hits);
            for (k, v) in pre_report{
                report.insert(colors_accession.get(&k).unwrap().to_string(), v);
            }
            if gene_search == false {
                super::generate_report(report, uniq_freqs, n_ref_kmers.to_owned(), cov);
            } else {
                super::generate_report_gene(report, num_kmers as usize);
            }
        } else {
//else we assume it is a fasta formatted file!
            println!("{}", file);
            eprintln!("Counting k-mers, this may take a while!");
            let vec_query = kmer_fa::read_fasta(file.to_owned());
            let unfiltered = kmer_fa::kmerize_vector(vec_query, k_size);
            let kmers_query = kmer_fa::clean_map(unfiltered, filter);
            let num_kmers = kmers_query.len() as f64;
            println!("{} k-mers in query", num_kmers);
            let mut report = HashMap::new();
            let mut uniq_freqs = HashMap::new();
            for (k, _) in &kmers_query {
                let mut kmer_slices = Vec::new();
                for i in 0..num_hash {
                    let bit_index = murmur_hash64a(k.as_bytes(), i as u64) % bloom_size as u64;
                    let bi = bit_index as usize;
                    if bigsi_map.contains_key(&bi) == false {
                        let count = report.entry(String::from("No hits!")).or_insert(0);
                        *count += 1;
                        break;
                    } else {
                        kmer_slices.push(bigsi_map.get(&bi).unwrap());
                    }
                }
                let mut first = BitVec::from_bytes(&kmer_slices[0].to_owned());
                for i in 1..num_hash {
                    let j = i as usize;
                    first.intersect(&BitVec::from_bytes(&kmer_slices[j]));
                }
                let mut hits = Vec::new();
                for i in 0..first.len() {
                    if first[i] == true {
                        hits.push(colors_accession.get(&i).unwrap());
                    }
                }
                for h in &hits {
                    let count = report.entry(h.to_string()).or_insert(0);
                    *count += 1;
                }
                if hits.len() == 1 {
                    let key = hits[0];
                    let value = *kmers_query.get(&k.to_string()).unwrap() as f64;
                    uniq_freqs
                        .entry(key.to_string())
                        .or_insert(Vec::new())
                        .push(value);
                }
            }
            if gene_search == false {
                super::generate_report(report, uniq_freqs, n_ref_kmers.to_owned(), cov);
            } else {
                super::generate_report_gene(report, num_kmers as usize);
            }
        }
    }
}
