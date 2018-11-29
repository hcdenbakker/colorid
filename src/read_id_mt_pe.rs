use bit_vec::BitVec;
use fasthash;
use flate2::read::MultiGzDecoder;
use fnv;
use kmer;
use probability::prelude::*;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use seq;
use std;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::sync::Arc;
use std::time::SystemTime;

pub fn false_prob_map(
    colors_accession: &fnv::FnvHashMap<usize, String>,
    ref_kmers_in: &fnv::FnvHashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
) -> fnv::FnvHashMap<usize, f64> {
    //first change a colors to accessions map to a accessions to colors map
    let mut accession_color = HashMap::new();
    for (color, accession) in colors_accession {
        accession_color.insert(accession, color);
    }
    let mut false_positive_p = fnv::FnvHashMap::default();
    for (key, value) in ref_kmers_in {
        let color = accession_color.get(key).unwrap();
        false_positive_p.insert(
            *color.to_owned(),
            false_prob(bloom_size as f64, num_hash as f64, *value as f64),
        );
    }
    false_positive_p
}

pub fn search_index(
    bigsi_map: &fnv::FnvHashMap<usize, BitVec>,
    map: &fnv::FnvHashMap<std::string::String, usize>,
    bloom_size: usize,
    num_hash: usize,
    no_hits_num: usize,
) -> fnv::FnvHashMap<usize, usize> {
    let mut report = fnv::FnvHashMap::default();
    let empty_bitvec = BitVec::from_elem(no_hits_num, false);
    for k in map.keys() {
        let mut kmer_slices = Vec::new();
        for i in 0..num_hash {
            let bit_index =
                fasthash::xx::hash64_with_seed(&k.as_bytes(), i as u64) % bloom_size as u64;
            let bi = bit_index as usize;
            if !bigsi_map.contains_key(&bi) || bigsi_map[&bi] == empty_bitvec {
                break;
            } else {
                kmer_slices.push(bigsi_map.get(&bi).unwrap());
            }
        }
        if kmer_slices.len() < num_hash {
            *report.entry(no_hits_num).or_insert(0) += 1;
            break;
        } else {
            let mut first = kmer_slices[0].to_owned();
            for i in 1..num_hash {
                let j = i as usize;
                first.intersect(&kmer_slices[j]);
            }
            let mut color = 0;
            for item in first {
                if item {
                    *report.entry(color).or_insert(0) += 1;
                }
                color += 1;
            }
        }
    }
    report
}

pub fn kmer_poll<'a>(
    report: &fnv::FnvHashMap<usize, usize>,
    _map: &fnv::FnvHashMap<std::string::String, usize>,
    colors_accession: &'a fnv::FnvHashMap<usize, String>,
    child_fp: &fnv::FnvHashMap<usize, f64>,
    no_hits_num: usize,
    fp_correct: f64,
) -> &'a str {
    let mut count_vec: Vec<_> = report.iter().collect();
    count_vec.sort_by(|a, b| b.1.cmp(a.1));
    let mut observations = 0;
    for (_k, v) in report {
        //if k == &no_hits_num {
        //    continue;
        //} else {
        observations += v;
        //}
    }
    if (count_vec[0].0 == &no_hits_num) && (count_vec.len() == 1) {
        //if observations == 0 {
        return "no_hits";
    };
    if count_vec[0].0 == &no_hits_num {
        let top_hit = count_vec[1].0;
        let p_false = child_fp.get(&top_hit).unwrap();
        let distribution = Binomial::new(observations, *p_false);
        let critical_value = observations as f64 * p_false;
        let mpf = distribution.mass(count_vec[1].1.to_owned());
        if ((count_vec[1].1.to_owned() as f64) < critical_value)
            || (((count_vec[1].1.to_owned() as f64) > critical_value) && (mpf >= fp_correct))
        {
            "reject"
        } else {
            &colors_accession[&top_hit]
        }
    } else {
        let top_hit = count_vec[0].0;
        let p_false = child_fp.get(&top_hit).unwrap();
        let distribution = Binomial::new(observations, *p_false);
        let critical_value = observations as f64 * p_false;
        let mpf = distribution.mass(count_vec[0].1.to_owned());
        if ((count_vec[0].1.to_owned() as f64) < critical_value)
            || (((count_vec[0].1.to_owned() as f64) > critical_value) && (mpf >= fp_correct))
        {
            "reject"
        } else {
            &colors_accession[&top_hit]
        }
    }
}

pub fn kmer_poll_mini<'a>(
    report: &fnv::FnvHashMap<usize, usize>,
    map: &fnv::FnvHashMap<std::string::String, usize>,
    colors_accession: &'a fnv::FnvHashMap<usize, String>,
    child_fp: &fnv::FnvHashMap<usize, f64>,
    no_hits_num: usize,
    fp_correct: f64,
) -> &'a str {
    let mut count_vec: Vec<_> = report.iter().collect();
    count_vec.sort_by(|a, b| b.1.cmp(a.1));
    let kmer_length = map.len();
    if (count_vec[0].0 == &no_hits_num) && (count_vec.len() == 1) {
        //if observations == 0 {
        return "no_hits";
    };
    if count_vec[0].0 == &no_hits_num {
        let top_hit = count_vec[1].0;
        let p_false = child_fp.get(&top_hit).unwrap();
        let distribution = Binomial::new(kmer_length, *p_false);
        let critical_value = kmer_length as f64 * p_false;
        let mpf = distribution.mass(count_vec[1].1.to_owned());
        if ((count_vec[1].1.to_owned() as f64) < critical_value)
            || (((count_vec[1].1.to_owned() as f64) > critical_value) && (mpf >= fp_correct))
        {
            "reject"
        } else {
            &colors_accession[&top_hit]
        }
    } else {
        let top_hit = count_vec[0].0;
        let p_false = child_fp.get(&top_hit).unwrap();
        let distribution = Binomial::new(kmer_length, *p_false);
        let critical_value = kmer_length as f64 * p_false;
        let mpf = distribution.mass(count_vec[0].1.to_owned());
        if ((count_vec[0].1.to_owned() as f64) < critical_value)
            || (((count_vec[0].1.to_owned() as f64) > critical_value) && (mpf >= fp_correct))
        {
            "reject"
        } else {
            &colors_accession[&top_hit]
        }
    }
}

//kmer poll classification plus raw count data as output
pub fn kmer_poll_plus<'a>(
    report: &fnv::FnvHashMap<usize, usize>,
    map: &fnv::FnvHashMap<std::string::String, usize>,
    colors_accession: &'a fnv::FnvHashMap<usize, String>,
    child_fp: &fnv::FnvHashMap<usize, f64>,
    no_hits_num: usize,
    fp_correct: f64,
) -> (&'a str, usize, usize, &'a str) {
    let mut count_vec: Vec<_> = report.iter().collect();
    count_vec.sort_by(|a, b| b.1.cmp(a.1));
    let kmer_length = map.len();
    if (count_vec[0].0 == &no_hits_num) && (count_vec.len() == 1) {
        return ("no_hits", 0 as usize, kmer_length, "accept");
    };
    if count_vec[0].0 == &no_hits_num {
        let top_hit = count_vec[1].0;
        let p_false = child_fp.get(&top_hit).unwrap();
        let distribution = Binomial::new(kmer_length, *p_false);
        let critical_value = kmer_length as f64 * p_false;
        let mpf = distribution.mass(count_vec[1].1.to_owned());
        if ((count_vec[1].1.to_owned() as f64) < critical_value)
            || (((count_vec[1].1.to_owned() as f64) > critical_value) && (mpf >= fp_correct))
        {
            (
                &colors_accession[&top_hit],
                count_vec[1].1.to_owned(),
                kmer_length,
                "reject",
            )
        } else {
            (
                &colors_accession[&top_hit],
                count_vec[1].1.to_owned(),
                kmer_length,
                "accept",
            )
        }
    } else {
        let top_hit = count_vec[0].0;
        let p_false = child_fp.get(&top_hit).unwrap();
        let distribution = Binomial::new(kmer_length, *p_false);
        let critical_value = kmer_length as f64 * p_false;
        let mpf = distribution.mass(count_vec[0].1.to_owned());
        if ((count_vec[0].1.to_owned() as f64) < critical_value)
            || (((count_vec[0].1.to_owned() as f64) > critical_value) && (mpf >= fp_correct))
        {
            (
                &colors_accession[&top_hit],
                count_vec[0].1.to_owned(),
                kmer_length,
                "reject",
            )
        } else {
            (
                &colors_accession[&top_hit],
                count_vec[0].1.to_owned(),
                kmer_length,
                "accept",
            )
        }
    }
}

pub struct Fasta {
    pub id: String,  //id including >
    pub seq: String, // sequence
}

impl Fasta {
    pub fn new() -> Fasta {
        Fasta {
            id: String::new(),
            seq: String::new(),
        }
    }
}

#[allow(unused_assignments)]
pub fn stream_fasta(
    filenames: Vec<&str>,
    bigsi_map: &fnv::FnvHashMap<usize, BitVec>, //has to be an Arc ?
    colors_accession: &fnv::FnvHashMap<usize, String>,
    ref_kmers_in: &fnv::FnvHashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
    k: usize,
    m: usize,
    t: usize, //threads
    d: usize, //downsample factor
    fp_correct: f64,
    b: usize,
    prefix: &str,
) /*-> std::collections::HashMap<std::string::String, usize>*/
{
    let f = File::open(filenames[0]).expect("file not found");
    let iter1 = io::BufReader::new(f).lines();
    //let mut tax_map = HashMap::new();
    let mut vec = Vec::with_capacity(b);
    //let mut master_vec = Vec::new();
    let no_hits_num: usize = colors_accession.len();
    let mut file =
        File::create(format!("{}_reads.txt", prefix)).expect("could not create outfile!");
    let mut sub_string = String::new();
    let mut fasta = Fasta::new();
    ThreadPoolBuilder::new()
        .num_threads(t)
        .build_global()
        .expect("Can't initialize ThreadPoolBuilder");
    let false_positive_p = false_prob_map(colors_accession, ref_kmers_in, bloom_size, num_hash);
    let my_bigsi: Arc<&fnv::FnvHashMap<usize, BitVec>> = Arc::new(bigsi_map);
    let false_positive_p_arc: Arc<fnv::FnvHashMap<usize, f64>> = Arc::new(false_positive_p);
    let search_time = SystemTime::now();
    let mut count = 0;
    let mut read_count = 0;
    for line in iter1 {
        let l = line.unwrap().to_string();
        if count == 0 {
            fasta.id = l;
        } else {
            if l.contains('>') {
                if sub_string.len() > 0 {
                    fasta.seq = sub_string.to_string();
                    vec.push(vec![fasta.id, fasta.seq]);
                    fasta.id = l;
                    sub_string.clear();
                }
            } else {
                sub_string.push_str(&l);
            }
        }
        count += 1;
        if vec.len() % b == 0 {
            let mut c: Vec<_> = vec![];
            c = vec
                .par_iter()
                .map(|r| {
                    let child_bigsi = my_bigsi.clone();
                    let child_fp = false_positive_p_arc.clone();
                    if r[1].len() < k {
                        (
                            r[0].to_owned(),
                            "too_short",
                            0 as usize,
                            0 as usize,
                            "accept",
                        )
                    } else {
                        let map = if m == 0 {
                            kmer::kmerize_vector_skip_n(vec![r[1].to_string()], k, d)
                        } else {
                            kmer::minimerize_vector_skip_n(vec![r[1].to_string()], k, m, d)
                        };
                        let report =
                            search_index(&child_bigsi, &map, bloom_size, num_hash, no_hits_num);
                        if report.is_empty() {
                            (r[0].to_owned(), "no_hits", 0 as usize, map.len(), "accept")
                        } else {
                            let classification = kmer_poll_plus(
                                &report,
                                &map,
                                &colors_accession,
                                &child_fp,
                                no_hits_num,
                                fp_correct,
                            );
                            (
                                r[0].to_owned(),
                                classification.0,
                                classification.1.to_owned(),
                                classification.2.to_owned(),
                                classification.3,
                            )
                        }
                    }
                }).collect();
            read_count += c.len();
            eprint!(" {} reads classified\r", read_count);
            for id in c {
                file.write_all(
                    format!("{}\t{}\t{}\t{}\t{}\n", id.0, id.1, id.2, id.3, id.4).as_bytes(),
                ).expect("could not write results!");
            }
            vec.clear();
        }
    }
    fasta.seq = sub_string.to_string();
    vec.push(vec![fasta.id, fasta.seq]);
    let mut c: Vec<_> = vec![];
    c = vec
        .par_iter()
        .map(|r| {
            let child_bigsi = my_bigsi.clone();
            let child_fp = false_positive_p_arc.clone();
            if r[1].len() < k {
                (
                    r[0].to_owned(),
                    "too_short",
                    0 as usize,
                    0 as usize,
                    "accept",
                )
            } else {
                let map = if m == 0 {
                    kmer::kmerize_vector_skip_n(vec![r[1].to_string()], k, d)
                } else {
                    kmer::minimerize_vector_skip_n(vec![r[1].to_string()], k, m, d)
                };
                let report = search_index(&child_bigsi, &map, bloom_size, num_hash, no_hits_num);
                if report.is_empty() {
                    (r[0].to_owned(), "no_hits", 0 as usize, map.len(), "accept")
                } else {
                    let classification = kmer_poll_plus(
                        &report,
                        &map,
                        &colors_accession,
                        &child_fp,
                        no_hits_num,
                        fp_correct,
                    );
                    (
                        r[0].to_owned(),
                        classification.0,
                        classification.1.to_owned(),
                        classification.2.to_owned(),
                        classification.3,
                    )
                }
            }
        }).collect();
    read_count += c.len();
    for id in c {
        file.write_all(format!("{}\t{}\t{}\t{}\t{}\n", id.0, id.1, id.2, id.3, id.4).as_bytes())
            .expect("could not write results!");
    }
    eprint!(" {} reads classified\r", read_count);
    match search_time.elapsed() {
        Ok(elapsed) => {
            eprintln!(
                "Classified {} reads in {} seconds",
                read_count,
                elapsed.as_secs()
            );
        }
        Err(e) => {
            // an error occurred!
            eprintln!("Error: {:?}", e);
        }
    }
}

pub fn false_prob(m: f64, k: f64, n: f64) -> f64 {
    let e = std::f64::consts::E;
    (1.0 - e.powf(-((k * (n + 0.5)) / (m - 1.0)))).powf(k)
}

#[allow(unused_assignments)]
pub fn per_read_stream_pe(
    filenames: Vec<&str>,
    bigsi_map: &fnv::FnvHashMap<usize, BitVec>, //has to be an Arc
    colors_accession: &fnv::FnvHashMap<usize, String>,
    ref_kmers_in: &fnv::FnvHashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
    k: usize,
    m: usize,
    t: usize, //threads
    d: usize, //downsample factor
    fp_correct: f64,
    b: usize,
    prefix: &str,
    qual_offset: u8,
) /*-> std::collections::HashMap<std::string::String, usize>*/
{
    let search_time = SystemTime::now();
    let mut vec = Vec::with_capacity(b);
    let mut line_count = 1;
    let no_hits_num: usize = colors_accession.len();
    let mut file =
        File::create(format!("{}_reads.txt", prefix)).expect("could not create outfile!");
    let batch = b * 4;
    let f = File::open(&filenames[0]).expect("file not found");
    let f2 = File::open(&filenames[1]).expect("file not found");
    let d1 = MultiGzDecoder::new(f);
    let d2 = MultiGzDecoder::new(f2);
    let iter1 = io::BufReader::new(d1).lines();
    let mut iter2 = io::BufReader::new(d2).lines();
    ThreadPoolBuilder::new()
        .num_threads(t)
        .build_global()
        .expect("Can't initialize ThreadPoolBuilder");
    let false_positive_p = false_prob_map(colors_accession, ref_kmers_in, bloom_size, num_hash);
    let my_bigsi: Arc<&fnv::FnvHashMap<usize, BitVec>> = Arc::new(bigsi_map);
    let false_positive_p_arc: Arc<fnv::FnvHashMap<usize, f64>> = Arc::new(false_positive_p);
    let mut read_count = 0;
    let _header = "".to_string();
    let mut fastq = seq::Fastq::new();
    for line in iter1 {
        let l = line.unwrap();
        let line2 = iter2.next();
        if line_count % 4 == 1 {
            fastq.id = l.to_string();
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
                    vec.push(vec![
                        fastq.id.to_owned(),
                        seq::qual_mask(fastq.seq1.to_owned(), l.to_owned(), qual_offset).to_owned(),
                        seq::qual_mask(fastq.seq2.to_owned(), l2.unwrap().to_owned(), qual_offset)
                            .to_owned(),
                    ]);
                }
                None => break,
            };
        }
        line_count += 1;
        if line_count % batch == 0 {
            let mut c: Vec<_> = vec![];
            c = vec
                .par_iter()
                .map(|r| {
                    let child_bigsi = my_bigsi.clone();
                    let child_fp = false_positive_p_arc.clone();
                    if (r[1].len() < k) || (r[2].len() < k) {
                        (
                            r[0].to_owned(),
                            "too_short",
                            0 as usize,
                            0 as usize,
                            "accept",
                        )
                    } else {
                        let map = if m == 0 {
                            kmer::kmerize_vector_skip_n(
                                vec![r[1].to_string(), r[2].to_string()],
                                k,
                                d,
                            )
                        } else {
                            kmer::minimerize_vector_skip_n(
                                vec![r[1].to_string(), r[2].to_string()],
                                k,
                                m,
                                d,
                            )
                        };
                        let report =
                            search_index(&child_bigsi, &map, bloom_size, num_hash, no_hits_num);
                        if report.is_empty() {
                            (r[0].to_owned(), "no_hits", 0 as usize, map.len(), "accept")
                        } else {
                            let classification = kmer_poll_plus(
                                &report,
                                &map,
                                &colors_accession,
                                &child_fp,
                                no_hits_num,
                                fp_correct,
                            );
                            (
                                r[0].to_owned(),
                                classification.0,
                                classification.1.to_owned(),
                                classification.2.to_owned(),
                                classification.3,
                            )
                        } //else
                    }
                }).collect();
            read_count += c.len();
            eprint!("{} read pairs classified\r", read_count);
            for id in c {
                file.write_all(
                    format!("{}\t{}\t{}\t{}\t{}\n", id.0, id.1, id.2, id.3, id.4).as_bytes(),
                ).expect("could not write results!");
            }
            vec.clear();
        }
    }
    //and the last block if there is any left
    let mut c: Vec<_> = vec![];
    c = vec
        .par_iter()
        .map(|r| {
            let child_bigsi = my_bigsi.clone();
            let child_fp = false_positive_p_arc.clone();
            if (r[1].len() < k) || (r[2].len() < k) {
                (
                    r[0].to_owned(),
                    "too_short",
                    0 as usize,
                    0 as usize,
                    "accept",
                )
            } else {
                let map = if m == 0 {
                    kmer::kmerize_vector_skip_n(vec![r[1].to_string(), r[2].to_string()], k, d)
                } else {
                    kmer::minimerize_vector_skip_n(
                        vec![r[1].to_string(), r[2].to_string()],
                        k,
                        m,
                        d,
                    )
                };
                let report = search_index(&child_bigsi, &map, bloom_size, num_hash, no_hits_num);
                if report.is_empty() {
                    (r[0].to_owned(), "no_hits", 0 as usize, map.len(), "accept")
                } else {
                    let classification = kmer_poll_plus(
                        &report,
                        &map,
                        &colors_accession,
                        &child_fp,
                        no_hits_num,
                        fp_correct,
                    );
                    (
                        r[0].to_owned(),
                        classification.0,
                        classification.1.to_owned(),
                        classification.2.to_owned(),
                        classification.3,
                    )
                } //else
            }
        }).collect();
    read_count += c.len();
    eprint!("{} read pairs classified\r", read_count);
    for id in c {
        file.write_all(format!("{}\t{}\t{}\t{}\t{}\n", id.0, id.1, id.2, id.3, id.4).as_bytes())
            .expect("could not write results!");
    }
    match search_time.elapsed() {
        Ok(elapsed) => {
            eprintln!(
                "Classified {} read pairs in {} seconds",
                read_count,
                elapsed.as_secs()
            );
        }
        Err(e) => {
            // an error occurred!
            eprintln!("Error: {:?}", e);
        }
    }
}

#[allow(unused_assignments)]
pub fn per_read_stream_se(
    filenames: Vec<&str>,
    bigsi_map: &fnv::FnvHashMap<usize, BitVec>, //has to be an Arc ?
    colors_accession: &fnv::FnvHashMap<usize, String>,
    ref_kmers_in: &fnv::FnvHashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
    k: usize,
    m: usize, //0 == no m, otherwise minimizer
    t: usize, //threads
    d: usize, //downsample factor
    fp_correct: f64,
    b: usize,
    prefix: &str,
    qual_offset: u8,
) /* -> std::collections::HashMap<std::string::String, usize>*/
{
    let search_time = SystemTime::now();
    let mut vec = Vec::with_capacity(b);
    let mut line_count = 1;
    let _header = "".to_string();
    let _sequence = "".to_string();
    let no_hits_num: usize = colors_accession.len();
    let mut read_count = 0;
    let mut file =
        File::create(format!("{}_reads.txt", prefix)).expect("could not create outfile!");
    let batch = b * 4;
    let _num_files = filenames.len();
    let f = File::open(&filenames[0]).expect("file not found");
    let d1 = MultiGzDecoder::new(f);
    let iter1 = io::BufReader::new(d1).lines();
    ThreadPoolBuilder::new()
        .num_threads(t)
        .build_global()
        .expect("Can't initialize ThreadPoolBuilder");
    let false_positive_p = false_prob_map(colors_accession, ref_kmers_in, bloom_size, num_hash);

    let my_bigsi: Arc<&fnv::FnvHashMap<usize, BitVec>> = Arc::new(bigsi_map);
    let false_positive_p_arc: Arc<fnv::FnvHashMap<usize, f64>> = Arc::new(false_positive_p);
    let mut fastq = seq::Fastq::new();
    for line in iter1 {
        let l = line.unwrap();
        if line_count % 4 == 1 {
            fastq.id = l.to_string();
        } else if line_count % 4 == 2 {
            fastq.seq1 = l.to_string();
        } else if line_count % 4 == 0 {
            vec.push(vec![
                fastq.id.to_owned(),
                fastq.seq1.to_owned(),
                l.to_owned(),
            ]);
        }
        line_count += 1;
        if line_count % batch == 0 {
            let mut c: Vec<_> = vec![];
            c = vec
                .par_iter()
                .map(|r| {
                    let child_bigsi = my_bigsi.clone();
                    let child_fp = false_positive_p_arc.clone();
                    if r[1].len() < k {
                        (
                            r[0].to_owned(),
                            "too_short",
                            0 as usize,
                            0 as usize,
                            "accept",
                        )
                    } else {
                        let map = if m == 0 {
                            kmer::kmerize_vector_skip_n(
                                vec![
                                    seq::qual_mask(r[1].to_string(), r[2].to_string(), qual_offset)
                                        .to_owned(),
                                ],
                                k,
                                d,
                            )
                        } else {
                            kmer::minimerize_vector_skip_n(
                                vec![
                                    seq::qual_mask(r[1].to_string(), r[2].to_string(), qual_offset)
                                        .to_owned(),
                                ],
                                k,
                                m,
                                d,
                            )
                        };
                        let report =
                            search_index(&child_bigsi, &map, bloom_size, num_hash, no_hits_num);
                        if report.is_empty() {
                            (r[0].to_owned(), "no_hits", 0 as usize, 0 as usize, "accept")
                        } else {
                            let classification = kmer_poll_plus(
                                &report,
                                &map,
                                &colors_accession,
                                &child_fp,
                                no_hits_num,
                                fp_correct,
                            );
                            (
                                r[0].to_owned(),
                                classification.0,
                                classification.1.to_owned(),
                                classification.2.to_owned(),
                                classification.3,
                            )
                        }
                    }
                }).collect();
            read_count += c.len();
            eprint!("{} read pairs classified\r", read_count);
            for id in c {
                file.write_all(
                    format!("{}\t{}\t{}\t{}\t{}\n", id.0, id.1, id.2, id.3, id.4).as_bytes(),
                ).expect("could not write results!");
            }
            vec.clear();
        }
    }
    //and the last block if there is any left
    let mut c: Vec<_> = vec![];
    c = vec
        .par_iter()
        .map(|r| {
            let child_bigsi = my_bigsi.clone();
            let child_fp = false_positive_p_arc.clone();
            if r[1].len() < k {
                (
                    r[0].to_owned(),
                    "too_short",
                    0 as usize,
                    0 as usize,
                    "accept",
                )
            } else {
                let map = if m == 0 {
                    kmer::kmerize_vector_skip_n(vec![r[1].to_string()], k, d)
                } else {
                    kmer::minimerize_vector_skip_n(vec![r[1].to_string()], k, m, d)
                };
                let report = search_index(&child_bigsi, &map, bloom_size, num_hash, no_hits_num);
                if report.is_empty() {
                    (r[0].to_owned(), "no_hits", 0 as usize, 0 as usize, "accept")
                } else {
                    let classification = kmer_poll_plus(
                        &report,
                        &map,
                        &colors_accession,
                        &child_fp,
                        no_hits_num,
                        fp_correct,
                    );
                    (
                        r[0].to_owned(),
                        classification.0,
                        classification.1.to_owned(),
                        classification.2.to_owned(),
                        classification.3,
                    )
                } //else
            }
        }).collect();
    read_count += c.len();
    eprint!("{} read pairs classified\r", read_count);
    for id in c {
        file.write_all(format!("{}\t{}\t{}\t{}\t{}\n", id.0, id.1, id.2, id.3, id.4).as_bytes())
            .expect("could not write results!");
    }
    match search_time.elapsed() {
        Ok(elapsed) => {
            eprintln!(
                "Classified {} reads in {} seconds",
                read_count,
                elapsed.as_secs()
            );
        }
        Err(e) => {
            // an error occurred!
            eprintln!("Error: {:?}", e);
        }
    }
}
