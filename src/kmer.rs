use flate2::read::MultiGzDecoder;
use fnv;
use seq;
use std;
use std::cmp;
use std::fs::File;
use std::io;
use std::io::prelude::*;

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

pub fn read_fasta_mf(filename: String) -> (Vec<String>, Vec<String>) {
    let mut f = File::open(filename).expect("file not found");
    let mut contents = String::new();
    f.read_to_string(&mut contents)
        .expect("something went wrong reading the file");
    let mut vec = Vec::new();
    let mut labels = Vec::new();
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
            labels.push(line[1..].to_owned());
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
    (labels, vec)
}

pub fn kmerize_vector(
    v: Vec<String>,
    k: usize,
    d: usize,
) -> fnv::FnvHashMap<std::string::String, usize> {
    let mut map = fnv::FnvHashMap::default();
    for l in v {
        let length_l = l.len();
        if length_l < k {
            continue;
        } else {
            let l_r = revcomp(&l);
            for i in 0..l.len() - k + 1 {
                if i % d == 0 {
                    if seq::has_no_n(l[i..i + k].as_bytes()) {
                        if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
                            let count = map
                                .entry(l[i..i + k].to_string().to_uppercase())
                                .or_insert(0);
                            *count += 1;
                        } else {
                            let count = map
                                .entry(
                                    l_r[length_l - (i + k)..length_l - i]
                                        .to_string()
                                        .to_uppercase(),
                                ).or_insert(0);
                            *count += 1;
                        }
                    }
                }
            }
        }
    }
    map
}

pub fn kmerize_vector_uppercase(
    v: Vec<String>,
    k: usize,
    d: usize,
) -> fnv::FnvHashMap<std::string::String, usize> {
    let mut map = fnv::FnvHashMap::default();
    for l in v {
        let length_l = l.len();
        let l_r = revcomp(&l);
        for i in 0..l.len() - k + 1 {
            if i % d == 0 {
                if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
                    let count = map.entry(l[i..i + k].to_string()).or_insert(0);
                    *count += 1;
                } else {
                    let count = map
                        .entry(l_r[length_l - (i + k)..length_l - i].to_string())
                        .or_insert(0);
                    *count += 1;
                }
            }
        }
    }
    map
}

pub fn kmerize_vector_skip_n(
    v: Vec<String>,
    k: usize,
    d: usize,
) -> fnv::FnvHashMap<std::string::String, usize> {
    let mut map = fnv::FnvHashMap::default();
    for l in v {
        let length_l = l.len();
        let l_r = revcomp(&l);
        for i in 0..l.len() - k + 1 {
            if i % d == 0 {
                if seq::has_no_n(l[i..i + k].as_bytes()) {
                    if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
                        let count = map.entry(l[i..i + k].to_string()).or_insert(0);
                        *count += 1;
                    } else {
                        let count = map
                            .entry(l_r[length_l - (i + k)..length_l - i].to_string())
                            .or_insert(0);
                        *count += 1;
                    }
                } else {
                    continue;
                }
            }
        }
    }
    map
}

pub fn kmerize_string(l: String, k: usize) -> fnv::FnvHashMap<std::string::String, usize> {
    let mut map = fnv::FnvHashMap::default();
    let length_l = l.len();
    let l_r = revcomp(&l);
    for i in 0..l.len() - k + 1 {
        if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
            let count = map
                .entry(l[i..i + k].to_string().to_uppercase())
                .or_insert(0);
            *count += 1;
        } else {
            let count = map
                .entry(
                    l_r[length_l - (i + k)..length_l - i]
                        .to_string()
                        .to_uppercase(),
                ).or_insert(0);
            *count += 1;
        }
    }
    map
}

pub fn minimerize_vector(
    v: Vec<String>,
    k: usize,
    m: usize,
    d: usize,
) -> fnv::FnvHashMap<std::string::String, usize> {
    let mut map = fnv::FnvHashMap::default();
    for l in v {
        let length_l = l.len();
        let l_r = revcomp(&l);
        for i in 0..l.len() - k + 1 {
            if i % d == 0 {
                if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
                    let min = find_minimizer(&l[i..i + k], m);
                    let count = map.entry(min).or_insert(0);
                    *count += 1;
                } else {
                    let min = find_minimizer(&l_r[length_l - (i + k)..length_l - i], m);
                    let count = map.entry(min).or_insert(0);
                    *count += 1;
                }
            }
        }
    }
    map
}

pub fn minimerize_vector_skip_n(
    v: Vec<String>,
    k: usize,
    m: usize,
    d: usize,
) -> fnv::FnvHashMap<std::string::String, usize> {
    let mut map = fnv::FnvHashMap::default();
    for l in v {
        let length_l = l.len();
        if length_l < k {
            continue;
        } else {
            let l_r = revcomp(&l);
            for i in 0..l.len() - k + 1 {
                if i % d == 0 {
                    if seq::has_no_n(l[i..i + k].as_bytes()) {
                        if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
                            let min = find_minimizer(&l[i..i + k], m);
                            let count = map.entry(min.to_uppercase()).or_insert(0);
                            *count += 1;
                        } else {
                            let min = find_minimizer(&l_r[length_l - (i + k)..length_l - i], m);
                            let count = map.entry(min.to_uppercase()).or_insert(0);
                            *count += 1;
                        }
                    } else {
                        continue;
                    }
                }
            }
        }
    }
    map
}

pub fn kmers_from_fq(filename: String, k: usize) -> fnv::FnvHashMap<std::string::String, usize> {
    let f = File::open(filename).expect("file not found");
    let mut map = fnv::FnvHashMap::default();
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
                        let count = map
                            .entry(l_r[length_l - (i + k)..length_l - i].to_string())
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

pub fn kmers_from_fq_qual(
    filename: String,
    k: usize,
    _d: usize,
    qual_offset: u8,
) -> fnv::FnvHashMap<std::string::String, usize> {
    let f = File::open(filename).expect("file not found");
    let gz = MultiGzDecoder::new(f);
    let iter = io::BufReader::new(gz).lines();
    let mut map = fnv::FnvHashMap::default();
    let mut line_count = 1;
    let mut fastq = seq::Fastq::new();
    let d = 1; //downsampler not used here so 1
    for line in iter {
        let l = line.unwrap();
        if line_count % 4 == 1 {
            fastq.id = l.to_string();
        } else if line_count % 4 == 2 {
            fastq.seq1 = l.to_owned();
        } else if line_count % 4 == 0 {
            fastq.qual1 = l.to_owned();
            let masked1 = seq::qual_mask(fastq.seq1.to_owned(), fastq.qual1, qual_offset);
            for l in vec![masked1] {
                let length_l = l.len();
                if length_l < k {
                    continue;
                } else {
                    let l_r = revcomp(&l);
                    for i in 0..l.len() - k + 1 {
                        if i % d == 0 {
                            if seq::has_no_n(l[i..i + k].as_bytes()) {
                                if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
                                    let count = map.entry(l[i..i + k].to_string()).or_insert(0);
                                    *count += 1;
                                } else {
                                    let count = map
                                        .entry(l_r[length_l - (i + k)..length_l - i].to_string())
                                        .or_insert(0);
                                    *count += 1;
                                }
                            }
                        }
                    }
                } //here
            }
        }
        line_count += 1;
    }
    map
}

pub fn kmers_from_fq_minimizer(
    filename: String,
    k: usize,
    m: usize,
) -> fnv::FnvHashMap<std::string::String, usize> {
    let f = File::open(filename).expect("file not found");
    let mut map = fnv::FnvHashMap::default();
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
                        let count = map.entry(find_minimizer(&l[i..i + k], m)).or_insert(0);
                        *count += 1;
                    } else {
                        let count = map
                            .entry(find_minimizer(&l_r[length_l - (i + k)..length_l - i], m))
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

pub fn kmers_fq_pe(filenames: Vec<&str>, k: usize) -> fnv::FnvHashMap<std::string::String, usize> {
    let mut map = fnv::FnvHashMap::default();
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
                            let count = map
                                .entry(l_r[length_l - (i + k)..length_l - i].to_string())
                                .or_insert(0);
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

pub fn kmers_fq_pe_qual(
    filenames: Vec<&str>,
    k: usize,
    d: usize,
    qual_offset: u8,
) -> fnv::FnvHashMap<std::string::String, usize> {
    let mut line_count = 1;
    let f = File::open(&filenames[0]).expect("file not found");
    let f2 = File::open(&filenames[1]).expect("file not found");
    let d1 = MultiGzDecoder::new(f);
    let d2 = MultiGzDecoder::new(f2);
    let iter1 = io::BufReader::new(d1).lines();
    let mut iter2 = io::BufReader::new(d2).lines();
    let mut map = fnv::FnvHashMap::default();
    let mut fastq = seq::Fastq::new();
    for line in iter1 {
        let l = line.unwrap();
        let line2 = iter2.next();
        if line_count % 4 == 1 {
            match line2 {
                Some(_h2) => {
                    fastq.id = l.to_string();
                }
                None => break,
            };
        } else if line_count % 4 == 2 {
            match line2 {
                Some(l2) => {
                    fastq.seq1 = l.to_owned();
                    fastq.seq2 = l2.unwrap().to_owned();
                }
                None => break,
            };
        } else if line_count % 4 == 0 {
            match line2 {
                Some(l2) => {
                    fastq.qual1 = l.to_owned();
                    fastq.qual2 = l2.unwrap().to_owned();
                    let masked1 = seq::qual_mask(fastq.seq1.to_owned(), fastq.qual1, qual_offset);
                    let masked2 = seq::qual_mask(fastq.seq2.to_owned(), fastq.qual2, qual_offset);
                    for l in vec![masked1, masked2] {
                        let length_l = l.len();
                        if length_l < k {
                            continue;
                        } else {
                            let l_r = revcomp(&l);
                            for i in 0..l.len() - k + 1 {
                                if i % d == 0 {
                                    if seq::has_no_n(l[i..i + k].as_bytes()) {
                                        if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
                                            let count =
                                                map.entry(l[i..i + k].to_string()).or_insert(0);
                                            *count += 1;
                                        } else {
                                            let count = map
                                                .entry(
                                                    l_r[length_l - (i + k)..length_l - i]
                                                        .to_string(),
                                                ).or_insert(0);
                                            *count += 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                None => break,
            };
        }
        line_count += 1;
    }
    map
}

pub fn kmers_fq_pe_minimizer(
    filenames: Vec<&str>,
    k: usize,
    m: usize,
) -> fnv::FnvHashMap<std::string::String, usize> {
    let mut map = fnv::FnvHashMap::default();
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
                            let count = map.entry(find_minimizer(&l[i..i + k], m)).or_insert(0);
                            *count += 1;
                        } else {
                            let count = map
                                .entry(find_minimizer(&l_r[length_l - (i + k)..length_l - i], m))
                                .or_insert(0);
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

pub fn kmers_fq_pe_minimizer_qual(
    filenames: Vec<&str>,
    k: usize,
    m: usize,
    d: usize,
    qual_offset: u8,
) -> fnv::FnvHashMap<std::string::String, usize> {
    let mut line_count = 1;
    let f = File::open(&filenames[0]).expect("file not found");
    let f2 = File::open(&filenames[1]).expect("file not found");
    let d1 = MultiGzDecoder::new(f);
    let d2 = MultiGzDecoder::new(f2);
    let iter1 = io::BufReader::new(d1).lines();
    let mut iter2 = io::BufReader::new(d2).lines();
    let mut map = fnv::FnvHashMap::default();
    let mut fastq = seq::Fastq::new();
    for line in iter1 {
        let l = line.unwrap();
        let line2 = iter2.next();
        if line_count % 4 == 1 {
            match line2 {
                Some(_h2) => {
                    fastq.id = l.to_string();
                }
                None => break,
            };
        } else if line_count % 4 == 2 {
            match line2 {
                Some(l2) => {
                    fastq.seq1 = l.to_owned();
                    fastq.seq2 = l2.unwrap().to_owned();
                }
                None => break,
            };
        } else if line_count % 4 == 0 {
            match line2 {
                Some(l2) => {
                    fastq.qual1 = l.to_owned();
                    fastq.qual2 = l2.unwrap().to_owned();
                    let masked1 = seq::qual_mask(fastq.seq1.to_owned(), fastq.qual1, qual_offset);
                    let masked2 = seq::qual_mask(fastq.seq2.to_owned(), fastq.qual2, qual_offset);
                    for l in vec![masked1, masked2] {
                        let length_l = l.len();
                        if length_l < k {
                            continue;
                        } else {
                            let l_r = revcomp(&l);
                            for i in 0..l.len() - k + 1 {
                                if i % d == 0 {
                                    if seq::has_no_n(l[i..i + k].as_bytes()) {
                                        if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
                                            let count = map
                                                .entry(find_minimizer(&l[i..i + k], m))
                                                .or_insert(0);
                                            *count += 1;
                                        } else {
                                            let count = map
                                                .entry(find_minimizer(
                                                    &l_r[length_l - (i + k)..length_l - i],
                                                    m,
                                                )).or_insert(0);
                                            *count += 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                None => break,
            };
        }
        line_count += 1;
    }
    map
}

pub fn kmers_from_fq_minimizer_qual(
    filename: String,
    k: usize,
    m: usize,
    _d: usize,
    qual_offset: u8,
) -> fnv::FnvHashMap<std::string::String, usize> {
    let f = File::open(filename).expect("file not found");
    let gz = MultiGzDecoder::new(f);
    let iter = io::BufReader::new(gz).lines();
    let mut map = fnv::FnvHashMap::default();
    let mut line_count = 1;
    let mut fastq = seq::Fastq::new();
    let d = 1; //downsampler not used here so 1
    for line in iter {
        let l = line.unwrap();
        if line_count % 4 == 1 {
            fastq.id = l.to_string();
        } else if line_count % 4 == 2 {
            fastq.seq1 = l.to_owned();
        } else if line_count % 4 == 0 {
            fastq.qual1 = l.to_owned();
            let masked1 = seq::qual_mask(fastq.seq1.to_owned(), fastq.qual1, qual_offset);
            for l in vec![masked1] {
                let length_l = l.len();
                if length_l < k {
                    continue;
                } else {
                    let l_r = revcomp(&l);
                    for i in 0..l.len() - k + 1 {
                        if i % d == 0 {
                            if seq::has_no_n(l[i..i + k].as_bytes()) {
                                if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
                                    let count =
                                        map.entry(find_minimizer(&l[i..i + k], m)).or_insert(0);
                                    *count += 1;
                                } else {
                                    let count = map
                                        .entry(find_minimizer(
                                            &l_r[length_l - (i + k)..length_l - i],
                                            m,
                                        )).or_insert(0);
                                    *count += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        line_count += 1;
    }
    map
}

pub fn clean_map(
    map: fnv::FnvHashMap<std::string::String, usize>,
    t: usize,
) -> fnv::FnvHashMap<std::string::String, usize> {
    let mut map_clean = fnv::FnvHashMap::default();
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

//auto cutoff inference from Zam Iqbal's Cortex
pub fn auto_cutoff(map: fnv::FnvHashMap<std::string::String, usize>) -> usize {
    let mut histo_map = fnv::FnvHashMap::default();
    for (_key, value) in map {
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

/* pseudominimizer based on murmurhash
pub fn find_minimizer(kmer: &str, m: usize) -> String {
    let mut minimizer = kmer[0..m].to_string();
    let mut value = murmur_hash64a(minimizer.as_bytes(), 0);
    for i in 1..kmer.len() - m + 1 {
        let alt_value = murmur_hash64a(kmer[i..i + m].as_bytes(), 0);
        if alt_value < value {
            minimizer = kmer[i..i + m].to_string();
            value = alt_value;
        }
    }
    minimizer
}*/

pub fn find_minimizer(kmer: &str, m: usize) -> String {
    let kmer_r = revcomp(&kmer);
    let length = kmer_r.len();
    let mut minimizers = Vec::new();
    for i in 1..kmer.len() - m + 1 {
        minimizers.push(kmer[i..i + m].to_string());
        minimizers.push(kmer_r[length - (i + m)..length - i].to_string());
    }
    minimizers.sort();
    minimizers[0].to_string()
}
