extern crate bincode;
extern crate bit_vec;
extern crate flate2;
extern crate kmer_fa;
extern crate murmurhash64;
extern crate serde;
#[macro_use]
extern crate serde_derive;
extern crate serde_json;
extern crate simple_bloom;

use std::collections::HashMap;
use std::io;
use std::io::prelude::*;
use std::fs::File;
use simple_bloom::BloomFilter;
use murmurhash64::murmur_hash64a;
use bit_vec::BitVec;
use serde::Serialize;
use serde_json::Serializer;
use std::io::Write;
use std::error::Error;
use flate2::read::GzDecoder;
use std::time::{Duration, SystemTime};

pub fn tab_to_map(filename: String) -> std::collections::HashMap<std::string::String, String> {
    let mut map = HashMap::new();
    let f = File::open(filename).expect("file not found");
    for line in io::BufReader::new(f).lines() {
        let l = line.unwrap();
        let v: Vec<&str> = l.split("\t").collect();
        map.insert(String::from(v[0]), String::from(v[1]));
    }
    map
}

pub fn build_bigsi2(
    map: std::collections::HashMap<std::string::String, String>,
    bloom_size: usize,
    num_hash: usize,
    k_size: usize,
) -> (
    std::collections::HashMap<usize, bit_vec::BitVec>,
    std::collections::HashMap<usize, String>,
    std::collections::HashMap<String, usize>,
) {
    //get a hash map with taxa and bit vector from bloom filter
    let mut bit_map = HashMap::new();
    let mut ref_kmer = HashMap::new();
    let mut accessions = Vec::new();
    for (accession, v) in &map {
        println!("Adding {} to BIGSI", accession);
        if v.ends_with("gz") {
            let unfiltered = kmer_fa::kmers_from_fq(v.to_owned(), k_size);
            let kmers = kmer_fa::clean_map(unfiltered, 1);
            let mut filter = BloomFilter::new(bloom_size as usize, num_hash as usize);
            for (kmer, _) in &kmers {
                filter.insert(&kmer);
            }
            bit_map.insert(accession, filter.bits);
        } else {
            let vec = kmer_fa::read_fasta(v.to_string());
            let kmers = kmer_fa::kmerize_vector(vec, k_size);
            ref_kmer.insert(accession.to_string(), kmers.len());
            let mut filter = BloomFilter::new(bloom_size as usize, num_hash as usize);
            for (kmer, _) in &kmers {
                filter.insert(&kmer);
            }
            bit_map.insert(accession, filter.bits);
        }
    }
    for (accession, _) in &map {
        accessions.push(accession);
    }
    //create hash table with colors for accessions
    let mut accession_colors = HashMap::new();
    let mut colors_accession = HashMap::new();
    for (c, s) in accessions.iter().enumerate() {
        accession_colors.insert(s.to_string(), c);
        colors_accession.insert(c, s.to_string());
    }
    let num_taxa = accessions.len();
    let mut bigsi_map = HashMap::new();
    for i in 0..bloom_size {
        let mut bitvec = bit_vec::BitVec::from_elem(num_taxa, false);
        for (t, s) in &bit_map {
            if s[i] == true {
                bitvec.set(*accession_colors.get(&t.to_string()).unwrap(), true);
            }
        }
        bigsi_map.insert(i, bitvec);
    }
    (bigsi_map, colors_accession, ref_kmer)
}

//there seems to be a ~18 s overhead to this function for larger BIGSIs
//add a hashmap with vectors containing k-mer coverages of uniquely placed k-mers to estimate
//coverage per taxon
pub fn search_bigsi(
    kmer_query: std::collections::HashMap<std::string::String, i32>,
    bigsi_map: std::collections::HashMap<usize, bit_vec::BitVec>,
    colors_accession: std::collections::HashMap<usize, String>,
    bloom_size: usize,
    num_hash: usize,
    k_size: usize,
) -> (
    std::collections::HashMap<String, usize>,
    std::collections::HashMap<String, Vec<f64>>,
) //hashmap with name and vector containing freqs taxon-specific kmers
{
    println!("Search! Collecting slices");
    let mut report = HashMap::new();
    let mut uniq_freqs = HashMap::new();
    for (k, _) in &kmer_query {
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
        let mut first = kmer_slices[0].to_owned();
        for i in 1..num_hash {
            let j = i as usize;
            first.intersect(&kmer_slices[j]);
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
            let value = *kmer_query.get(&k.to_string()).unwrap() as f64;
            uniq_freqs
                .entry(key.to_string())
                .or_insert(Vec::new())
                .push(value);
        }
    }
    (report, uniq_freqs)
}
// test with search_bigsi function show this will  introduce a ~13 second time overhead as opposed
// to writing the search out
pub fn batch_search(
    files: Vec<&str>,
    bigsi_map: std::collections::HashMap<usize, bit_vec::BitVec>,
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
            println!("Counting k-mers, this may take a while!");
            let unfiltered = kmer_fa::kmers_from_fq(file.to_owned(), k_size);
            let kmers_query = kmer_fa::clean_map(unfiltered, filter);
            let num_kmers = kmers_query.len() as f64;
            println!("{} k-mers in query", num_kmers);
            let bigsi_search = SystemTime::now();
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
                let mut first = kmer_slices[0].to_owned();
                for i in 1..num_hash {
                    let j = i as usize;
                    first.intersect(&kmer_slices[j]);
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
            match bigsi_search.elapsed() {
                Ok(elapsed) => {
                    println!("Search: {} sec", elapsed.as_secs());
                }
                Err(e) => {
                    // an error occurred!
                    println!("Error: {:?}", e);
                }
            }
            if gene_search == false {
                generate_report(report, uniq_freqs, n_ref_kmers.to_owned(), cov);
            } else {
                generate_report_gene(report, num_kmers as usize);
            }
        } else {
            //else we assume it is a fasta formatted file!
            println!("{}", file);
            println!("Counting k-mers, this may take a while!");
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
                let mut first = kmer_slices[0].to_owned();
                for i in 1..num_hash {
                    let j = i as usize;
                    first.intersect(&kmer_slices[j]);
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
                generate_report(report, uniq_freqs, n_ref_kmers.to_owned(), cov);
            } else {
                generate_report_gene(report, num_kmers as usize);
            }
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct BigsyMap {
    pub bloom_size: usize,
    pub num_hash: usize,
    pub k_size: usize,
    pub colors: HashMap<usize, String>,
    pub map: HashMap<usize, Vec<u8>>,
    pub n_ref_kmers: HashMap<String, usize>,
}

pub fn save_bigsi(
    bigsi: std::collections::HashMap<usize, bit_vec::BitVec>,
    colors_accession: std::collections::HashMap<usize, String>,
    n_ref_kmers_in: std::collections::HashMap<String, usize>,
    bloom_size_in: usize,
    num_hash_in: usize,
    k_size_in: usize,
    path: &str,
)
/*-> Result<(), Box<Error>>*/
{
    let mut bigsi_map = HashMap::new();
    for (k, v) in bigsi {
        bigsi_map.insert(k, v.to_bytes());
    }
    let mappy = BigsyMap {
        map: bigsi_map,
        colors: colors_accession,
        n_ref_kmers: n_ref_kmers_in,
        bloom_size: bloom_size_in,
        num_hash: num_hash_in,
        k_size: k_size_in,
    };
    let serialized = serde_json::to_string(&mappy).unwrap();
    let mut writer = File::create(path).unwrap();
    writer.write_all(serialized.as_bytes());
}

pub fn read_bigsi(
    path: &str,
) -> (
    std::collections::HashMap<usize, bit_vec::BitVec>,
    std::collections::HashMap<usize, String>,
    std::collections::HashMap<String, usize>,
    usize,
    usize,
    usize,
) {
    let mut file = File::open(path).unwrap();
    let mut contents = String::new();
    file.read_to_string(&mut contents);
    let deserialized: BigsyMap = serde_json::from_str(&contents).unwrap();
    let mut bigsi_map = HashMap::new();
    //create a full hashmap with bitvectors (first naive implementation, takes care of 'empty'
    //bitslices...may not be necessary...
    //for i in 0..deserialized.bloom_size{
    //    bigsi_map.insert(i, bit_vec::BitVec::from_elem(deserialized.colors.len(),false));
    //}
    for (key, vector) in deserialized.map {
        bigsi_map.insert(key, BitVec::from_bytes(&vector));
    }
    (
        bigsi_map,
        deserialized.colors,
        deserialized.n_ref_kmers,
        deserialized.bloom_size,
        deserialized.num_hash,
        deserialized.k_size,
    )
}

//https://codereview.stackexchange.com/questions/173338/calculate-mean-median-and-mode-in-rust
pub fn mode(numbers: &Vec<f64>) -> usize {
    let mut occurrences = HashMap::new();

    for value in numbers {
        *occurrences.entry(*value as usize).or_insert(0) += 1;
    }

    occurrences
        .into_iter()
        .max_by_key(|&(_, count)| count)
        .map(|(val, _)| val)
        .expect("Cannot compute the mode of zero numbers")
}

pub fn generate_report(
    report: std::collections::HashMap<String, usize>,
    uniq_freqs: std::collections::HashMap<String, Vec<f64>>,
    n_ref_kmers: std::collections::HashMap<String, usize>,
    cov: f64,
) {
    for (k, v) in report {
        let frequencies = uniq_freqs.get(&k.to_string());
        let mut mean: f64 = 0.0;
        let mut modus: usize = 0;
        let mut specific_kmers: usize;
        match frequencies {
            Some(_x) => {
                mean = frequencies.unwrap().iter().fold(0.0, |a, &b| a + b)
                    / frequencies.unwrap().len() as f64;
                modus = mode(frequencies.unwrap());
                specific_kmers = frequencies.unwrap().len();
            }
            None => {
                mean = 0.0;
                modus = 0;
                specific_kmers = 0;
            }
        }
        let n_kmers = n_ref_kmers.get(&k.to_string());
        match n_kmers {
            Some(_x) => {
                let genome_cov = v as f64 / *n_kmers.unwrap() as f64;
                if genome_cov > cov {
                    println!(
                        "{}: {:.2} {:.2} {} {}",
                        k, genome_cov, mean, modus, specific_kmers
                    );
                }
            }
            None => continue,
        }
    }
}

pub fn generate_report_gene(
    report: std::collections::HashMap<String, usize>,
    gene_kmer_size: usize,
) {
    for (k, v) in report {
        let gene_match = v as f64 / gene_kmer_size as f64;
        if gene_match > 0.35 {
            println!("{}: {:.2}", k, gene_match);
        }
    }
}
