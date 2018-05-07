extern crate flate2;
extern crate murmurhash64;
extern crate probability;  
extern crate kmer_fa;
extern crate bit_vec;
extern crate rayon;

use rayon::prelude::*;
use std;
use std::sync::Arc;
use std::io;
use std::io::prelude::*;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use probability::prelude::*;
use std::collections::HashMap;
use murmurhash64::murmur_hash64a;
use bit_vec::BitVec;

pub fn nuc_reads_from_fq(filename: &str) -> Vec<String>{
   let mut f = File::open(filename).expect("file not found");
   let mut vec = Vec::new();
   let mut line_count = 1;
   let d = MultiGzDecoder::new(f);
   for line in io::BufReader::new(d).lines() {
       let l =  line.unwrap();
       if line_count%4 == 2 {
        vec.push(l);
        }

        line_count += 1;
    }

    vec
}

pub fn false_prob(m: f64, k: f64, n: f64) -> f64 {
    let e = std::f64::consts::E;
    let prob = (1.0 - e.powf(-((k * (n + 0.5)) / (m - 1.0)))).powf(k);
    prob
}


pub fn per_read_search(
    filename: String,
    bigsi_map: std::collections::HashMap<usize, Vec<u8>>,//has to be an Arc
    colors_accession: std::collections::HashMap<usize, String>,//has to be an Arc
    ref_kmers_in: std::collections::HashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
    k: usize,
    t: usize,
) -> std::collections::HashMap<std::string::String, usize> {
    rayon::ThreadPoolBuilder::new().num_threads(t).build_global().unwrap();
    let mut false_positive_p = HashMap::new();
    for (key, value) in ref_kmers_in {
        false_positive_p.insert(
            key,
            false_prob(bloom_size as f64, num_hash as f64, value as f64),
        );
    }
    let my_bigsi: Arc<std::collections::HashMap<usize, Vec<u8>>> = Arc::new(bigsi_map);
    //let my_colors: Arc<std::collections::HashMap<usize, String>> = Arc::new(colors_accession);
    let false_positive_p_Arc: Arc<std::collections::HashMap<std::string::String, f64>> = Arc::new(false_positive_p.clone());
    let reads = nuc_reads_from_fq(&filename);
    let mut tax_map = HashMap::new();
    let mut c: Vec<_> = vec![];
    c = reads.par_iter().map(|l|
            {
            let child_bigsi = my_bigsi.clone();
            let child_fp = false_positive_p_Arc.clone();
            let mut map = HashMap::new();
            let mut report = HashMap::new();
            let l_r = kmer_fa::revcomp(&l);
            let length_l = l.len();
            if length_l < k {
                return "too_short";
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
                for (k, _) in &map {
                    let mut kmer_slices = Vec::new();
                    for i in 0..num_hash {
                        let bit_index = murmur_hash64a(k.as_bytes(), i as u64) % bloom_size as u64;
                                                let bi = bit_index as usize;
                        if (child_bigsi.contains_key(&bi) == false) || (child_bigsi[&bi] == empty_bitvec){
                            *report.entry("No hits!").or_insert(0) += 1;
                            break;
                        } else {
                            kmer_slices.push(child_bigsi.get(&bi).unwrap());
                        }
                    }
                    let mut first = BitVec::from_bytes(&kmer_slices[0].to_owned());
                    for i in 1..num_hash {
                        let j = i as usize;
                        first.intersect(&BitVec::from_bytes(&kmer_slices[j]));
                    }
                    let mut color = 0;
                    for i in first{
                        if i ==true{
                            *report.entry(colors_accession.get(&color).unwrap()).or_insert(0) += 1;
                        }
                        color += 1;
                    }
                }
                let kmer_length = length_l - k + 1;
                let mut count_vec: Vec<_> = report.iter().collect();
                count_vec.sort_by(|a, b| b.1.cmp(a.1));
                if count_vec.len() == 0 {
                    return "no_hits";
                } else {
                    let top_hit = count_vec[0].0;
                    let p_false = child_fp.get(*top_hit).unwrap();
                    let distribution = Binomial::new(kmer_length, *p_false);
                    let critical_value = kmer_length as f64 * p_false;
                    let mpf = distribution.mass(count_vec[0].1.to_owned());
                    if (count_vec[0].1.to_owned() as f64) < critical_value {
                        return "no_hits";
                    } else if ((count_vec[0].1.to_owned() as f64) > critical_value)
                        && (mpf >= 0.001)
                    {
                        return "no_hits";
                    } else {
                        return top_hit;
                    }
                }
            }
    }).collect();
    eprint!("Classified {} reads\r", c.len());
    eprint!("\n");
    for id in c{
        let count = tax_map.entry(id.to_string()).or_insert(0);
            *count += 1;
    }
    tax_map
}
