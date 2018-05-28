extern crate flate2;

use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::cmp;

//add a fasta to kmer bloom filter

pub fn read_fasta(filename: String) -> Vec<String> {
    let mut f = File::open(filename).expect("file not found");
    let mut contents = String::new();
    f.read_to_string(&mut contents)
        .expect("something went wrong reading the file");
    let mut vec = Vec::new();
    let mut sub_string: String = "".to_owned();
    let mut vec_raw = Vec::new();
    for line in contents.lines() {
        let l = line.to_string();
        vec_raw.push(l);
    }
    let length_raw_vec = vec_raw.len();
    let mut count_line = 0;
    for line in vec_raw {
        count_line += 1;
        if line.contains('>') {
            let l = sub_string.to_string();
            if l.len() > 0 {
                vec.push(l);
            }
            sub_string.clear();
        } else if count_line == length_raw_vec {
            let l = line.to_string();
            sub_string.push_str(&l);
            let l = sub_string.to_string();
            if l.len() > 0 {
                vec.push(l);
            }
        } else {
            let l = line.to_string();
            sub_string.push_str(&l);
        }
    }
    vec
}

pub fn kmerize_vector(
    v: Vec<String>,
    k: usize,
) -> std::collections::HashMap<std::string::String, usize> {
    let mut map = HashMap::new();
    for l in v {
        let length_l = l.len();
        let l_r = revcomp(&l);
        for i in 0..l.len() - k + 1 {
            if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
                let count = map.entry(l[i..i + k].to_string().to_uppercase())
                    .or_insert(0);
                *count += 1;
            } else {
                let count = map.entry(
                    l_r[length_l - (i + k)..length_l - i]
                        .to_string()
                        .to_uppercase(),
                ).or_insert(0);
                *count += 1;
            }
        }
    }
    map
}

pub fn kmers_from_fq(
    filename: String,
    k: usize,
) -> std::collections::HashMap<std::string::String, usize> {
    let mut f = File::open(filename).expect("file not found");
    let mut map = HashMap::new();
    let mut line_count = 1;
    let d = MultiGzDecoder::new(f);
    for line in io::BufReader::new(d).lines() {
        let l = line.unwrap();
        let length_l = l.len();
        if line_count % 4 == 2 {
            if length_l < k {
                continue;
            } else {
                let l_r = revcomp(&l);
                for i in 0..l.len() - k + 1 {
                    if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
                        let count = map.entry(l[i..i + k].to_string()).or_insert(0);
                        *count += 1;
                    } else {
                        let count = map.entry(l_r[length_l - (i + k)..length_l - i].to_string())
                            .or_insert(0);
                        *count += 1;
                    }
                }
            }
        }
        line_count += 1;
    }

    map
}

pub fn kmers_fq_pe(
    filenames: Vec<&str>,
    k: usize,
) -> std::collections::HashMap<std::string::String, usize> {
    let mut map = HashMap::new();
    for filename in filenames {
        let mut f = File::open(filename).expect("file not found");
        let mut line_count = 1;
        let d = MultiGzDecoder::new(f);
        for line in io::BufReader::new(d).lines() {
            let l = line.unwrap();
            let length_l = l.len();
            if line_count % 4 == 2 {
                if length_l < k {
                    continue;
                } else {
                    let l_r = revcomp(&l);
                    for i in 0..l.len() - k + 1 {
                        if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
                            let count = map.entry(l[i..i + k].to_string()).or_insert(0);
                            *count += 1;
                        } else {
                            let count = map.entry(
                                l_r[length_l - (i + k)..length_l - i].to_string(),
                            ).or_insert(0);
                            *count += 1;
                        }
                    }
                }
            }
            line_count += 1;
        }
    }
    map
}

//auto cutoff inference from Zam Iqbal's Cortex
pub fn auto_cutoff(map: std::collections::HashMap<std::string::String, usize>) -> usize {
    let mut histo_map = HashMap::new();
    for (key, value) in map {
        *histo_map.entry(value).or_insert(0) += 1;
    }
    let mut count_vec: Vec<_> = histo_map.iter().collect();
    count_vec.sort_by(|&(a, _), &(b, _)| a.cmp(&b));
    let mut coverages: Vec<usize> = Vec::with_capacity(count_vec.len());
    for v in count_vec {
        coverages.push(*v.1);
    }
    //first pseudo-derivative
    let mut d1 = Vec::new();
    for i in 1..coverages.len() - 1 {
        d1.push(coverages[i] as f64 / coverages[i + 1] as f64);
    }
    //second pseudo-derivative
    let mut d2 = Vec::new();
    for i in 0..d1.len() - 1 {
        d2.push(d1[i] / d1[i + 1]);
    }
    let mut first_pos_d1 = 0;
    let mut first_pos_d2 = 0;
    let threshold: f64 = 1.0;
    for (i, p) in d1.iter().enumerate() {
        if p < &threshold {
            first_pos_d1 = i + 1;
            break;
        }
    }
    for (i, p) in d2.iter().enumerate() {
        if p < &threshold {
            first_pos_d2 = i + 1;
            break;
        }
    }
    //estimate coverage (mean), exclude singleton k-mers
    let mut bigsum = 0;
    for (i, p) in coverages[1..].iter().enumerate() {
        bigsum += i * p;
    }
    let num_kmers: usize = coverages[1..].iter().sum();
    let mean: f64 = bigsum as f64 / num_kmers as f64;
    if (first_pos_d1 > 0) && ((first_pos_d1 as f64) < (mean * 0.75)) {
        first_pos_d1
    } else if first_pos_d2 > 0 {
        first_pos_d2
    } else {
        cmp::max(1, (mean / 2.0).ceil() as usize)
    }
}

pub fn clean_map(
    map: std::collections::HashMap<std::string::String, usize>,
    t: usize,
) -> std::collections::HashMap<std::string::String, usize> {
    let mut map_clean = HashMap::new();
    for (key, value) in map {
        if value > t {
            map_clean.insert(key, value);
        }
    }
    map_clean
}

pub fn revcomp(dna: &str) -> String {
    let mut rc_dna: String = String::with_capacity(dna.len());
    for c in dna.chars().rev() {
        rc_dna.push(switch_base(c))
    }
    rc_dna
}

fn switch_base(c: char) -> char {
    match c {
        'a' => 't',
        'c' => 'g',
        't' => 'a',
        'g' => 'c',
        'u' => 'a',
        'n' => 'n',
        'A' => 'T',
        'C' => 'G',
        'T' => 'A',
        'G' => 'C',
        'U' => 'A',
        'N' => 'N',
        _ => 'N',
    }
}
