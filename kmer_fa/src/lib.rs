extern crate flate2;

use std::collections::HashMap;
use std::io;
use std::io::prelude::*;
use flate2::read::GzDecoder;
use std::fs::File;

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
        if line.contains('>') || count_line == length_raw_vec {
            let l = sub_string.to_string();
            if l.len() > 0 {
                vec.push(l);
            }
            sub_string.clear();
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
) -> std::collections::HashMap<std::string::String, i32> {
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
) -> std::collections::HashMap<std::string::String, i32> {
    let mut f = File::open(filename).expect("file not found");
    let mut map = HashMap::new();
    let mut line_count = 1;
    let d = GzDecoder::new(f);
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



pub fn clean_map(
    map: std::collections::HashMap<std::string::String, i32>,
    t: i32,
) -> std::collections::HashMap<std::string::String, i32> {
    let mut map_clean = HashMap::new();
    for (key, value) in map {
        if value > t {
            map_clean.insert(key, value);
        }
    }
    map_clean
}

fn revcomp(dna: &str) -> String {
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
