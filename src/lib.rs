extern crate bincode;
extern crate bit_vec;
extern crate flate2;
extern crate kmer_fa;
extern crate murmurhash64;
extern crate probability;
extern crate rayon;
extern crate serde;
#[macro_use]
extern crate serde_derive;
extern crate serde_json;
extern crate simple_bloom;

use bincode::{deserialize, serialize, Infinite};
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

//adjust for pe implementation
pub fn tab_to_map_old(filename: String) -> std::collections::HashMap<std::string::String, String> {
    let mut map = HashMap::new();
    let f = File::open(filename).expect("reference file not found");
    for line in io::BufReader::new(f).lines() {
        let l = line.unwrap();
        let v: Vec<&str> = l.split('\t').collect();
        map.insert(String::from(v[0]), String::from(v[1]));
    }
    map
}

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

pub fn build_bigsi2(
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
            let unfiltered = kmer_fa::kmers_fq_pe(vec![&v[0], &v[1]], k_size);
            let cutoff = kmer_fa::auto_cutoff(unfiltered.to_owned());
            let kmers = kmer_fa::clean_map(unfiltered, cutoff);
            let mut filter = BloomFilter::new(bloom_size as usize, num_hash as usize);
            for kmer in kmers.keys() {
                filter.insert(&kmer);
            }
            bit_map.insert(accession, filter.bits);
        } else {
            if v[0].ends_with("gz") {
                let unfiltered = kmer_fa::kmers_from_fq(v[0].to_owned(), k_size);
                let kmers = kmer_fa::clean_map(unfiltered, 1);
                let mut filter = BloomFilter::new(bloom_size as usize, num_hash as usize);
                for kmer in kmers.keys() {
                    filter.insert(&kmer);
                }
                bit_map.insert(accession, filter.bits);
            } else {
                let vec = kmer_fa::read_fasta(v[0].to_string());
                let kmers = kmer_fa::kmerize_vector(vec, k_size);
                ref_kmer.insert(accession.to_string(), kmers.len());
                let mut filter = BloomFilter::new(bloom_size as usize, num_hash as usize);
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

//there seems to be a ~18 s overhead to this function for larger BIGSIs
//add a hashmap with vectors containing k-mer coverages of uniquely placed k-mers to estimate
//coverage per taxon
pub fn search_bigsi(
    kmer_query: &std::collections::HashMap<std::string::String, i32>,
    bigsi_map: &std::collections::HashMap<usize, Vec<u8>>,
    colors_accession: &std::collections::HashMap<usize, String>,
    bloom_size: usize,
    num_hash: usize,
) -> (
    std::collections::HashMap<String, usize>,
    std::collections::HashMap<String, Vec<f64>>,
    std::collections::HashMap<String, Vec<f64>>,
) {
    eprintln!("Search! Collecting slices");
    let mut report = HashMap::new();
    let mut uniq_freqs = HashMap::new();
    let mut multi_freqs = HashMap::new();
    for kmer in kmer_query.keys() {
        let mut kmer_slices = Vec::new();
        for i in 0..num_hash {
            let bit_index = murmur_hash64a(kmer.as_bytes(), i as u64) % bloom_size as u64;
            let bi = bit_index as usize;
            if !bigsi_map.contains_key(&bi) {
                let count = report.entry(String::from("No hits!")).or_insert(0);
                *count += 1;
                break;
            } else {
                kmer_slices.push(&bigsi_map[&bi]);
            }
        }
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
            let value = f64::from(kmer_query[&kmer.to_string()]);
            multi_freqs
                .entry(h.to_string())
                .or_insert_with(Vec::new)
                .push(value);
        }
        if hits.len() == 1 {
            let key = hits[0];
            let value = f64::from(kmer_query[&kmer.to_string()]);
            uniq_freqs
                .entry(key.to_string())
                .or_insert_with(Vec::new)
                .push(value);
        }
    }
    (report, uniq_freqs, multi_freqs)
}

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
    filter: isize,
    cov: f64,
    gene_search: bool,
) {
    for file in files {
        if file.ends_with("gz") {
            eprintln!("{}", file);
            eprintln!("Counting k-mers, this may take a while!");
            let unfiltered = kmer_fa::kmers_from_fq(file.to_owned(), k_size);
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
                generate_report(file, report, &uniq_freqs, &n_ref_kmers, cov);
            } else {
                generate_report_gene(file, report, num_kmers as usize);
            }
        } else {
            //else we assume it is a fasta formatted file!
            eprintln!("{}", file);
            eprintln!("Counting k-mers, this may take a while!");
            let vec_query = kmer_fa::read_fasta(file.to_owned());
            let unfiltered = kmer_fa::kmerize_vector(vec_query, k_size);
            let kmers_query =
                if filter < 0{
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
                generate_report(file, report, &uniq_freqs, &n_ref_kmers, cov);
            } else {
                generate_report_gene(file, report, num_kmers as usize);
            }
        }
    }
}

pub mod build_mt;

#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct BigsyMap {
    pub bloom_size: usize,
    pub num_hash: usize,
    pub k_size: usize,
    pub colors: HashMap<usize, String>,
    pub map: HashMap<usize, Vec<u8>>,
    pub n_ref_kmers: HashMap<String, usize>,
}

pub mod read_id_mt_v3;

pub mod read_id_mt_pe;

pub mod perfect_search;

pub mod batch_search_pe;

pub fn save_bigsi(
    bigsi_map: std::collections::HashMap<usize, Vec<u8>>,
    colors_accession: std::collections::HashMap<usize, String>,
    n_ref_kmers_in: std::collections::HashMap<String, usize>,
    bloom_size_in: usize,
    num_hash_in: usize,
    k_size_in: usize,
    path: &str,
) {
    let mappy = BigsyMap {
        map: bigsi_map,
        colors: colors_accession,
        n_ref_kmers: n_ref_kmers_in,
        bloom_size: bloom_size_in,
        num_hash: num_hash_in,
        k_size: k_size_in,
    };
    let serialized: Vec<u8> = serialize(&mappy, Infinite).unwrap();
    let mut writer = File::create(path).unwrap();
    writer
        .write_all(&serialized)
        .expect("problems preparing serialized data for writing");
}

pub fn save_bigsi_gz(
    bigsi_map: std::collections::HashMap<usize, Vec<u8>>,
    colors_accession: std::collections::HashMap<usize, String>,
    n_ref_kmers_in: std::collections::HashMap<String, usize>,
    bloom_size_in: usize,
    num_hash_in: usize,
    k_size_in: usize,
    path: &str,
) {
    let mappy = BigsyMap {
        map: bigsi_map,
        colors: colors_accession,
        n_ref_kmers: n_ref_kmers_in,
        bloom_size: bloom_size_in,
        num_hash: num_hash_in,
        k_size: k_size_in,
    };
    let serialized: Vec<u8> = serialize(&mappy, Infinite).unwrap();
    let f = File::create(path).unwrap();
    let mut gz = GzEncoder::new(f, Compression::new(9));
    gz.write_all(&serialized)
        .expect("could not write serialized index to file");
    gz.finish().expect("did not finish writing file properly");
}

pub fn read_bigsi(
    path: &str,
) -> (
    std::collections::HashMap<usize, Vec<u8>>,
    std::collections::HashMap<usize, String>,
    std::collections::HashMap<String, usize>,
    usize,
    usize,
    usize,
) {
    let mut reader = BufReader::new(File::open(path).expect("Can't open index!"));
    let mut buffer = Vec::new();
    reader.read_to_end(&mut buffer).expect("Can't read content");
    let deserialized: BigsyMap = deserialize(&buffer[..]).expect("cant deserialize");
    (
        deserialized.map,
        deserialized.colors,
        deserialized.n_ref_kmers,
        deserialized.bloom_size,
        deserialized.num_hash,
        deserialized.k_size,
    )
}

pub fn read_bigsi_gz(
    path: &str,
) -> (
    std::collections::HashMap<usize, Vec<u8>>,
    std::collections::HashMap<usize, String>,
    std::collections::HashMap<String, usize>,
    usize,
    usize,
    usize,
) {
    let file = File::open(path).expect("file not found");
    let mut gz = GzDecoder::new(file);
    let mut contents = Vec::new();
    gz.read_to_end(&mut contents)
        .expect("could not read contents file");
    let deserialized: BigsyMap = deserialize(&contents[..]).expect("cant deserialize");
    /*let mut bigsi_map = IndexMap::with_capacity(deserialized.map.len());
    for (key, vector) in deserialized.map {
        bigsi_map.insert(key, BitVec::from_bytes(&vector));
    }*/
    (
        deserialized.map,
        deserialized.colors,
        deserialized.n_ref_kmers,
        deserialized.bloom_size,
        deserialized.num_hash,
        deserialized.k_size,
    )
}
//https://codereview.stackexchange.com/questions/173338/calculate-mean-median-and-mode-in-rust
pub fn mode(numbers: &[f64]) -> usize {
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
    query: &str,
    report: std::collections::HashMap<String, usize>,
    uniq_freqs: &std::collections::HashMap<String, Vec<f64>>,
    n_ref_kmers: &std::collections::HashMap<String, usize>,
    cov: f64,
) {
    for (k, v) in report {
        let frequencies = uniq_freqs.get(&k.to_string());
        let mut mean: f64;
        let mut modus: usize;
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
                        "{}\t{}\t{:.2}\t{:.2}\t{}\t{}",
                        query, k, genome_cov, mean, modus, specific_kmers
                    );
                }
            }
            None => continue,
        }
    }
}

pub fn generate_report_plus(
    report: std::collections::HashMap<String, usize>,
    uniq_freqs: &std::collections::HashMap<String, Vec<f64>>,
    multi_freqs: &std::collections::HashMap<String, Vec<f64>>,
    n_ref_kmers: &std::collections::HashMap<String, usize>,
    cov: f64,
) {
    for (k, v) in report {
        let frequencies = uniq_freqs.get(&k.to_string());
        let mut mean: f64;
        let mut modus: usize;
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
        let multi_frequencies = multi_freqs.get(&k.to_string());
        let mut multi_mean: f64;
        let mut multi_modus: usize;
        let mut multi_kmers: usize;
        match multi_frequencies {
            Some(_x) => {
                multi_mean = multi_frequencies.unwrap().iter().fold(0.0, |a, &b| a + b)
                    / multi_frequencies.unwrap().len() as f64;
                multi_modus = mode(multi_frequencies.unwrap());
                multi_kmers = multi_frequencies.unwrap().len();
            }
            None => {
                multi_mean = 0.0;
                multi_modus = 0;
                multi_kmers = 0;
            }
        }
        let n_kmers = n_ref_kmers.get(&k.to_string());
        match n_kmers {
            Some(_x) => {
                let genome_cov = v as f64 / *n_kmers.unwrap() as f64;
                if genome_cov > cov {
                    println!(
                        "{}: {:.2} {:.2} {} {:.2} {} {}",
                        k, genome_cov, multi_mean, multi_modus, mean, modus, specific_kmers
                    );
                }
            }
            None => continue,
        }
    }
}

pub fn generate_report_gene(
    query: &str,
    report: std::collections::HashMap<String, usize>,
    gene_kmer_size: usize,
) {
    for (k, v) in report {
        let gene_match = v as f64 / gene_kmer_size as f64;
        if gene_match > 0.35 {
            println!("{}\t{}\t{}\t{:.3}", query, k, gene_kmer_size, gene_match);
        }
    }
}

pub fn false_prob(m: f64, k: f64, n: f64) -> f64 {
    let e = std::f64::consts::E;
    (1.0 - e.powf(-((k * (n + 0.5)) / (m - 1.0)))).powf(k)
}
