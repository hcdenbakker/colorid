use bit_vec::BitVec;
use flate2::read::MultiGzDecoder;
use kmer;
use murmurhash64::murmur_hash64a;
use probability::prelude::*;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::sync::Arc;
use std::time::SystemTime;

pub fn nuc_reads_from_fq(filename: &str) -> Vec<Vec<String>> {
    let f = File::open(filename).expect("file not found");
    let mut vec = Vec::new();
    let mut line_count = 1;
    let d = MultiGzDecoder::new(f);
    for line in io::BufReader::new(d).lines() {
        let l = line.unwrap();
        if line_count % 4 == 2 {
            vec.push(vec![l]);
        }

        line_count += 1;
    }

    vec
}

pub fn read_fasta(filename: String) -> Vec<Vec<String>> {
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
                vec.push(vec![l]);
            }
            sub_string.clear();
        } else if count_line == length_raw_vec {
            let l = line.to_string();
            sub_string.push_str(&l);
            let l = sub_string.to_string();
            if l.len() > 0 {
                vec.push(vec![l]);
            }
        } else {
            let l = line.to_string();
            sub_string.push_str(&l);
        }
    }
    vec
}

pub fn nuc_reads_fq_pe(filename1: &str, filename2: &str) -> Vec<Vec<String>> {
    let f1 = File::open(filename1).expect("file 1 not found");
    let f2 = File::open(filename2).expect("file 2 not found");
    let mut vec1 = Vec::new();
    let mut vec2 = Vec::new();
    //read forward reads
    let mut line_count1 = 1;
    let d1 = MultiGzDecoder::new(f1);
    for line in io::BufReader::new(d1).lines() {
        let l = line.unwrap();
        if line_count1 % 4 == 2 {
            vec1.push(l);
        }

        line_count1 += 1;
    }
    //read reverse reads
    let mut line_count2 = 1;
    let d2 = MultiGzDecoder::new(f2);
    for line in io::BufReader::new(d2).lines() {
        let l = line.unwrap();
        if line_count2 % 4 == 2 {
            vec2.push(l);
        }

        line_count2 += 1;
    }
    //check if forward and reverse have the same length
    if vec1.len() != vec2.len() {
        panic!("Forward and reverse file contain different number of reads!")
    }
    let mut vec_combined = Vec::new();
    for (i, r) in vec1.iter().enumerate() {
        vec_combined.push(vec![r.to_owned(), vec2[i].to_owned()])
    }
    vec_combined
}

pub fn false_prob(m: f64, k: f64, n: f64) -> f64 {
    let e = std::f64::consts::E;
    (1.0 - e.powf(-((k * (n + 0.5)) / (m - 1.0)))).powf(k)
}

pub fn per_read_search(
    filenames: Vec<&str>,
    bigsi_map: std::collections::HashMap<usize, Vec<u8>>, //has to be an Arc
    colors_accession: &std::collections::HashMap<usize, String>,
    ref_kmers_in: &std::collections::HashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
    k: usize,
    t: usize,
    d: usize, //downsample factor
    fp_correct: f64,
) -> std::collections::HashMap<std::string::String, usize> {
    ThreadPoolBuilder::new()
        .num_threads(t)
        .build_global()
        .unwrap();
    let mut false_positive_p = HashMap::new();
    for (key, value) in ref_kmers_in {
        false_positive_p.insert(
            key.to_owned(),
            false_prob(bloom_size as f64, num_hash as f64, *value as f64),
        );
    }
    let my_bigsi: Arc<std::collections::HashMap<usize, Vec<u8>>> = Arc::new(bigsi_map);
    let false_positive_p_arc: Arc<std::collections::HashMap<std::string::String, f64>> =
        Arc::new(false_positive_p.clone());
    let reads: Vec<Vec<String>> =
        //paired end
        if filenames.len() == 2{
            nuc_reads_fq_pe(&filenames[0], &filenames[1])
        }else{
        if filenames[0].ends_with("gz") {
        nuc_reads_from_fq(&filenames[0])
    } else {
        read_fasta(filenames[0].to_owned())
    }
    };
    let search_time = SystemTime::now();
    let mut tax_map = HashMap::new();
    let c: Vec<_>;
    c = reads
        .par_iter()
        .map(|r| {
            let child_bigsi = my_bigsi.clone();
            let child_fp = false_positive_p_arc.clone();
            let mut map = HashMap::new();
            let mut report = HashMap::new();
            if (r.len() == 1) && (r[0].len() < k) {
                "too_short"
            } else if (r[0].len() < k) && (r[1].len() < k) {
                "too_short"
            } else {
                for l in r {
                    let l_r = kmer::revcomp(&l);
                    let length_l = l.len();
                    if length_l < k {
                        break;
                    } else {
                        for i in 0..l.len() - k + 1 {
                            if i % d == 0 {
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
                }
                let empty_bitvec = BitVec::from_elem(child_fp.len(), false).to_bytes();
                for k in map.keys() {
                    let mut kmer_slices = Vec::new();
                    for i in 0..num_hash {
                        let bit_index = murmur_hash64a(k.as_bytes(), i as u64) % bloom_size as u64;
                        let bi = bit_index as usize;
                        if !child_bigsi.contains_key(&bi) || (child_bigsi[&bi] == empty_bitvec) {
                            break;
                        } else {
                            kmer_slices.push(child_bigsi.get(&bi).unwrap());
                        }
                    }
                    if kmer_slices.len() < num_hash {
                        //eprintln!("short");
                        *report.entry("no_hits").or_insert(0) += 1;
                        break;
                    } else {
                        let mut first = BitVec::from_bytes(&kmer_slices[0].to_owned());
                        for i in 1..num_hash {
                            let j = i as usize;
                            first.intersect(&BitVec::from_bytes(&kmer_slices[j]));
                        }
                        let mut color = 0;
                        for item in first {
                            if item {
                                *report.entry(&colors_accession[&color]).or_insert(0) += 1;
                            }
                            color += 1;
                        }
                    }
                }
                if report.is_empty() {
                    "no_hits"
                } else {
                    let mut count_vec: Vec<_> = report.iter().collect();
                    count_vec.sort_by(|a, b| b.1.cmp(a.1));
                    let kmer_length = map.len();
                    let top_hit = count_vec[0].0;
                    if *top_hit == "no_hits" {
                        return "no_hits";
                    }
                    let p_false = child_fp.get(*top_hit).unwrap();
                    let distribution = Binomial::new(kmer_length, *p_false);
                    let critical_value = kmer_length as f64 * p_false;
                    let mpf = distribution.mass(count_vec[0].1.to_owned());
                    if ((count_vec[0].1.to_owned() as f64) < critical_value)
                        || (((count_vec[0].1.to_owned() as f64) > critical_value)
                            && (mpf >= fp_correct))
                    {
                        "no_hits"
                    } else {
                        top_hit
                    }
                }
            }
        })
        .collect();
    match search_time.elapsed() {
        Ok(elapsed) => {
            eprintln!(
                "Classified {} reads in {} seconds",
                c.len(),
                elapsed.as_secs()
            );
        }
        Err(e) => {
            // an error occurred!
            eprintln!("Error: {:?}", e);
        }
    }
    for id in c {
        let count = tax_map.entry(id.to_string()).or_insert(0);
        *count += 1;
    }
    tax_map
}

pub fn per_read_search_minimizer(
    filenames: Vec<&str>,
    bigsi_map: std::collections::HashMap<usize, Vec<u8>>, //has to be an Arc
    colors_accession: &std::collections::HashMap<usize, String>,
    ref_kmers_in: &std::collections::HashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
    k: usize,
    m: usize, //minimizer size
    t: usize, //threads
    d: usize, //downsample factor
    fp_correct: f64,
) -> std::collections::HashMap<std::string::String, usize> {
    ThreadPoolBuilder::new()
        .num_threads(t)
        .build_global()
        .unwrap();
    let mut false_positive_p = HashMap::new();
    for (key, value) in ref_kmers_in {
        false_positive_p.insert(
            key.to_owned(),
            false_prob(bloom_size as f64, num_hash as f64, *value as f64),
        );
    }
    let my_bigsi: Arc<std::collections::HashMap<usize, Vec<u8>>> = Arc::new(bigsi_map);
    let false_positive_p_arc: Arc<std::collections::HashMap<std::string::String, f64>> =
        Arc::new(false_positive_p.clone());
    let reads: Vec<Vec<String>> =
        //paired end
        if filenames.len() == 2{
            nuc_reads_fq_pe(&filenames[0], &filenames[1])
        }else{
        if filenames[0].ends_with("gz") {
        nuc_reads_from_fq(&filenames[0])
    } else {
        read_fasta(filenames[0].to_owned())
    }
    };
    let search_time = SystemTime::now();
    let mut tax_map = HashMap::new();
    let c: Vec<_>;
    c = reads
        .par_iter()
        .map(|r| {
            let child_bigsi = my_bigsi.clone();
            let child_fp = false_positive_p_arc.clone();
            let mut map = HashMap::new();
            let mut report = HashMap::new();
            if (r.len() == 1) && (r[0].len() < k) {
                "too_short"
            } else if (r[0].len() < k) && (r[1].len() < k) {
                "too_short"
            } else {
                for l in r {
                    let l_r = kmer::revcomp(&l);
                    let length_l = l.len();
                    if length_l < k {
                        break;
                    } else {
                        for i in 0..l.len() - k + 1 {
                            if i % d == 0 {
                                if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
                                    let count = map.entry(kmer::find_minimizer(&l[i..i + k], m))
                                        .or_insert(0);
                                    *count += 1;
                                } else {
                                    let count = map.entry(
                                        kmer::find_minimizer(
                                            &l_r[length_l - (i + k)..length_l - i],
                                            m,
                                        ).to_string(),
                                    ).or_insert(0);
                                    *count += 1;
                                }
                            }
                        }
                    }
                }
                let empty_bitvec = BitVec::from_elem(child_fp.len(), false).to_bytes();
                for k in map.keys() {
                    let mut kmer_slices = Vec::new();
                    for i in 0..num_hash {
                        let bit_index = murmur_hash64a(k.as_bytes(), i as u64) % bloom_size as u64;
                        let bi = bit_index as usize;
                        if !child_bigsi.contains_key(&bi) || (child_bigsi[&bi] == empty_bitvec) {
                            break;
                        } else {
                            kmer_slices.push(child_bigsi.get(&bi).unwrap());
                        }
                    }
                    if kmer_slices.len() < num_hash {
                        //eprintln!("short");
                        *report.entry("no_hits").or_insert(0) += 1;
                        break;
                    } else {
                        let mut first = BitVec::from_bytes(&kmer_slices[0].to_owned());
                        for i in 1..num_hash {
                            let j = i as usize;
                            first.intersect(&BitVec::from_bytes(&kmer_slices[j]));
                        }
                        let mut color = 0;
                        for item in first {
                            if item {
                                *report.entry(&colors_accession[&color]).or_insert(0) += 1;
                            }
                            color += 1;
                        }
                    }
                }
                if report.is_empty() {
                    "no_hits"
                } else {
                    let mut count_vec: Vec<_> = report.iter().collect();
                    count_vec.sort_by(|a, b| b.1.cmp(a.1));
                    let kmer_length = map.len();
                    let top_hit = count_vec[0].0;
                    if *top_hit == "no_hits" {
                        return "no_hits";
                    }
                    let p_false = child_fp.get(*top_hit).unwrap();
                    let distribution = Binomial::new(kmer_length, *p_false);
                    let critical_value = kmer_length as f64 * p_false;
                    let mpf = distribution.mass(count_vec[0].1.to_owned());
                    if ((count_vec[0].1.to_owned() as f64) < critical_value)
                        || (((count_vec[0].1.to_owned() as f64) > critical_value)
                            && (mpf >= fp_correct))
                    {
                        "no_hits"
                    } else {
                        top_hit
                    }
                }
            }
        })
        .collect();
    match search_time.elapsed() {
        Ok(elapsed) => {
            eprintln!(
                "Classified {} reads in {} seconds",
                c.len(),
                elapsed.as_secs()
            );
        }
        Err(e) => {
            // an error occurred!
            eprintln!("Error: {:?}", e);
        }
    }
    for id in c {
        let count = tax_map.entry(id.to_string()).or_insert(0);
        *count += 1;
    }
    tax_map
}
