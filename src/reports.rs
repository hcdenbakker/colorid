use std;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::prelude::*;
use fnv;

pub fn generate_report(
    query: &str,
    report: std::collections::HashMap<String, usize>,
    uniq_freqs: &std::collections::HashMap<String, Vec<f64>>,
    n_ref_kmers: &fnv::FnvHashMap<String, usize>,
    num_kmers: usize,
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
                        "{}\t{}\t{}\t{:.2}\t{:.2}\t{}\t{}",
                        query, num_kmers, k, genome_cov, mean, modus, specific_kmers
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
    cov: f64,
) {
    for (k, v) in report {
        let gene_match = v as f64 / gene_kmer_size as f64;
        if gene_match >= cov {
            println!("{}\t{}\t{}\t{:.3}", query, k, gene_kmer_size, gene_match);
        }
    }
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

pub fn read_counts(file_name: String, prefix: &str) {
    let f = File::open(file_name).expect("file not found");
    let iter1 = io::BufReader::new(f).lines();
    let mut count_map = fnv::FnvHashMap::default();
    for line in iter1 {
        let l = line.unwrap().to_string();
        let v: Vec<&str> = l.split('\t').collect();
        let count = count_map.entry(v[1].to_owned()).or_insert(0);
        *count += 1;
    }
    let mut count_file =
        File::create(format!("{}_counts.txt", prefix)).expect("could not create outfile!");
    for (key, value) in &count_map {
        count_file
            .write_all(format!("{}\t {}\n", key, value).as_bytes())
            .expect("could not write count results!");
    }
}
