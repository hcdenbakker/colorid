extern crate bit_vec;
extern crate flate2;
extern crate kmer_fa;
extern crate murmurhash64;
extern crate probability;
extern crate rayon;

use bit_vec::BitVec;
use flate2::read::MultiGzDecoder;
use murmurhash64::murmur_hash64a;
use probability::prelude::*;
use rayon::prelude::*;
use std;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::sync::Arc;
use std::time::SystemTime;

pub fn nuc_reads_from_fq(filename: &str) -> Vec<String> {
    let f = File::open(filename).expect("file not found");
    let mut vec = Vec::new();
    let mut line_count = 1;
    let d = MultiGzDecoder::new(f);
    for line in io::BufReader::new(d).lines() {
        let l = line.unwrap();
        if line_count % 4 == 2 {
            vec.push(l);
        }

        line_count += 1;
    }

    vec
}

pub fn false_prob(m: f64, k: f64, n: f64) -> f64 {
    let e = std::f64::consts::E;
    (1.0 - e.powf(-((k * (n + 0.5)) / (m - 1.0)))).powf(k)
}

pub fn per_read_search(
    filename: String,
    bigsi_map: std::collections::HashMap<usize, Vec<u8>>, //has to be an Arc
    colors_accession: &std::collections::HashMap<usize, String>,
    ref_kmers_in: &std::collections::HashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
    k: usize,
    t: usize,
    fp_correct: f64,
) -> std::collections::HashMap<std::string::String, usize> {
    rayon::ThreadPoolBuilder::new()
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
    let reads: Vec<String> = if filename.ends_with("gz") {
        nuc_reads_from_fq(&filename)
    } else {
        kmer_fa::read_fasta(filename)
    };
    let search_time = SystemTime::now();
    let mut tax_map = HashMap::new();
    let c: Vec<_>;
    c = reads
        .par_iter()
        .map(|l| {
            let child_bigsi = my_bigsi.clone();
            let child_fp = false_positive_p_arc.clone();
            let mut map = HashMap::new();
            let mut report = HashMap::new();
            let l_r = kmer_fa::revcomp(&l);
            let length_l = l.len();
            if length_l < k {
                "too_short"
            } else {
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
                let empty_bitvec = bit_vec::BitVec::from_elem(child_fp.len(), false).to_bytes();
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
                    let kmer_length = length_l - k + 1;
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
