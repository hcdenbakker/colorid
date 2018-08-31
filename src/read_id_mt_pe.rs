use bit_vec::BitVec;
use flate2::read::MultiGzDecoder;
use kmer;
use kmer::minimerize_vector;
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

pub fn false_prob_map(
    colors_accession: &std::collections::HashMap<usize, String>,
    ref_kmers_in: &std::collections::HashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
) -> std::collections::HashMap<usize, f64> {
    //first change a colors to accessions map to a accessions to colors map
    let mut accession_color = HashMap::new();
    for (color, accession) in colors_accession {
        accession_color.insert(accession, color);
    }
    let mut false_positive_p = HashMap::new();
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
    bigsi_map: &std::collections::HashMap<usize, Vec<u8>>,
    map: &std::collections::HashMap<std::string::String, usize>,
    bloom_size: usize,
    num_hash: usize,
    no_hits_num: usize,
) -> std::collections::HashMap<usize, usize> {
    let mut report = HashMap::new();
    let empty_bitvec = BitVec::from_elem(no_hits_num, false).to_bytes();
    for k in map.keys() {
        let mut kmer_slices = Vec::new();
        for i in 0..num_hash {
            let bit_index = murmur_hash64a(k.as_bytes(), i as u64) % bloom_size as u64;
            let bi = bit_index as usize;
            if !bigsi_map.contains_key(&bi) || bigsi_map[&bi] == empty_bitvec {
                break;
            } else {
                kmer_slices.push(bigsi_map.get(&bi).unwrap());
            }
        }
        if kmer_slices.len() < num_hash {
            //eprintln!("short");
            *report.entry(no_hits_num).or_insert(0) += 1;
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
                    *report.entry(color).or_insert(0) += 1;
                }
                color += 1;
            }
        }
    }
    report
}

pub fn kmer_pol<'a>(
    report: &std::collections::HashMap<usize, usize>,
    map: &std::collections::HashMap<std::string::String, usize>,
    colors_accession: &'a std::collections::HashMap<usize, String>,
    child_fp: &std::collections::HashMap<usize, f64>,
    no_hits_num: usize,
    fp_correct: f64,
) -> &'a str {
    let mut count_vec: Vec<_> = report.iter().collect();
    count_vec.sort_by(|a, b| b.1.cmp(a.1));
    let kmer_length = map.len();
    let top_hit = count_vec[0].0;
    if *top_hit == no_hits_num {
        return "no_hits";
    }
    let p_false = child_fp.get(&top_hit).unwrap();
    let distribution = Binomial::new(kmer_length, *p_false);
    let critical_value = kmer_length as f64 * p_false;
    let mpf = distribution.mass(count_vec[0].1.to_owned());
    if ((count_vec[0].1.to_owned() as f64) < critical_value)
        || (((count_vec[0].1.to_owned() as f64) > critical_value) && (mpf >= fp_correct))
    {
        "no_hits"
    } else {
        &colors_accession[&top_hit]
    }
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
    let search_time = SystemTime::now();
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
    match search_time.elapsed() {
        Ok(elapsed) => {
            eprintln!(
                "Read {} read pairs in {} seconds",
                vec_combined.len(),
                elapsed.as_secs()
            );
        }
        Err(e) => {
            // an error occurred!
            eprintln!("Error: {:?}", e);
        }
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
    m: usize,
    t: usize,
    d: usize, //downsample factor
    fp_correct: f64,
) -> std::collections::HashMap<std::string::String, usize> {
    ThreadPoolBuilder::new()
        .num_threads(t)
        .build_global()
        .unwrap();
    let no_hits_num: usize = colors_accession.len();
    //let mut false_positive_p = HashMap::new();
    let false_positive_p = false_prob_map(colors_accession, ref_kmers_in, bloom_size, num_hash);
    let my_bigsi: Arc<std::collections::HashMap<usize, Vec<u8>>> = Arc::new(bigsi_map);
    let false_positive_p_arc: Arc<std::collections::HashMap<usize, f64>> =
        Arc::new(false_positive_p);
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
            if (r.len() == 1) && (r[0].len() < k) {
                "too_short"
            } else if (r[0].len() < k) && (r[1].len() < k) {
                "too_short"
            } else {
                let map = if m == 0 {
                    kmer::kmerize_vector(r.to_vec(), k, d)
                } else {
                    kmer::minimerize_vector(r.to_vec(), k, m, d)
                };
                let report = search_index(&child_bigsi, &map, bloom_size, num_hash, no_hits_num);
                if report.is_empty() {
                    "no_hits"
                } else {
                    kmer_pol(
                        &report,
                        &map,
                        &colors_accession,
                        &child_fp,
                        no_hits_num,
                        fp_correct,
                    )
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

pub fn per_read_stream_pe(
    filenames: Vec<&str>,
    bigsi_map: &std::collections::HashMap<usize, Vec<u8>>, //has to be an Arc ?
    colors_accession: &std::collections::HashMap<usize, String>,
    ref_kmers_in: &std::collections::HashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
    k: usize,
    m: usize,
    t: usize, //threads
    d: usize, //downsample factor
    fp_correct: f64,
    b: usize,
    prefix: &str,
) -> std::collections::HashMap<std::string::String, usize> {
    let search_time = SystemTime::now();
    let mut tax_map = HashMap::new();
    let mut vec = Vec::with_capacity(b);
    let mut master_vec = Vec::new();
    let mut line_count = 1;
    let no_hits_num: usize = colors_accession.len();
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
        .unwrap();
    let false_positive_p = false_prob_map(colors_accession, ref_kmers_in, bloom_size, num_hash);
    let my_bigsi: Arc<&std::collections::HashMap<usize, Vec<u8>>> = Arc::new(bigsi_map);
    let false_positive_p_arc: Arc<std::collections::HashMap<usize, f64>> =
        Arc::new(false_positive_p);
    let mut header = "".to_string();
    for line in iter1 {
        let l = line.unwrap();
        let line2 = iter2.next();
        if line_count % 4 == 1 {
            header = l.to_string();
        } else if line_count % 4 == 2 {
            match line2 {
                Some(l2) => {
                    vec.push(vec![
                        header.to_owned(),
                        l.to_owned(),
                        l2.unwrap().to_owned(),
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
                    if (r.len() == 2) && (r[1].len() < k) {
                        (r[0].to_owned(), "too_short")
                    } else if (r[1].len() < k) && (r[2].len() < k) {
                        (r[0].to_owned(), "too_short")
                    } else {
                        let map = if m == 0 {
                            kmer::kmerize_vector(vec![r[1].to_string(), r[2].to_string()], k, d)
                        } else {
                            kmer::minimerize_vector(
                                vec![r[1].to_string(), r[2].to_string()],
                                k,
                                m,
                                d,
                            )
                        };
                        let report =
                            search_index(&child_bigsi, &map, bloom_size, num_hash, no_hits_num);
                        if report.is_empty() {
                            (r[0].to_owned(), "no_hits")
                        } else {
                            (
                                r[0].to_owned(),
                                kmer_pol(
                                    &report,
                                    &map,
                                    &colors_accession,
                                    &child_fp,
                                    no_hits_num,
                                    fp_correct,
                                ),
                            )
                        }
                    }
                })
                .collect();
            master_vec.append(&mut c);
            eprint!("{} read pairs classified\r", master_vec.len());
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
            if (r[1].len() < k) && (r[4].len() < k) {
                (r[0].to_owned(), "too_short")
            } else {
                let map = if m == 0 {
                    kmer::kmerize_vector(vec![r[1].to_string(), r[2].to_string()], k, d)
                } else {
                    kmer::minimerize_vector(vec![r[1].to_string(), r[2].to_string()], k, m, d)
                };
                let report = search_index(&child_bigsi, &map, bloom_size, num_hash, no_hits_num);
                if report.is_empty() {
                    (r[0].to_owned(), "no_hits")
                } else {
                    (
                        r[0].to_owned(),
                        kmer_pol(
                            &report,
                            &map,
                            &colors_accession,
                            &child_fp,
                            no_hits_num,
                            fp_correct,
                        ),
                    )
                }
            }
        })
        .collect();
    master_vec.append(&mut c);
    eprintln!("{} read pairs classified\n", master_vec.len());
    match search_time.elapsed() {
        Ok(elapsed) => {
            eprintln!(
                "Classified {} read pairs in {} seconds",
                master_vec.len(),
                elapsed.as_secs()
            );
        }
        Err(e) => {
            // an error occurred!
            eprintln!("Error: {:?}", e);
        }
    }
    let mut file =
        File::create(format!("{}_reads.txt", prefix)).expect("could not create outfile!");
    for id in master_vec {
        file.write_all(format!("{}\t {}\n", id.0, id.1).as_bytes())
            .expect("could not write results!");
        let count = tax_map.entry(id.1.to_string()).or_insert(0);
        *count += 1;
    }
    let mut count_file =
        File::create(format!("{}_counts.txt", prefix)).expect("could not create outfile!");
    for (k, v) in &tax_map {
        count_file
            .write_all(format!("{}\t {}\n", k, v).as_bytes())
            .expect("could not write count results!");
    }
    tax_map
}

pub fn per_read_stream_se(
    filenames: Vec<&str>,
    bigsi_map: &std::collections::HashMap<usize, Vec<u8>>, //has to be an Arc ?
    colors_accession: &std::collections::HashMap<usize, String>,
    ref_kmers_in: &std::collections::HashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
    k: usize,
    m: usize, //0 == no m, otherwise minimizer
    t: usize, //threads
    d: usize, //downsample factor
    fp_correct: f64,
    b: usize,
    prefix: &str,
) -> std::collections::HashMap<std::string::String, usize> {
    let search_time = SystemTime::now();
    let mut tax_map = HashMap::new();
    let mut vec = Vec::with_capacity(b);
    let mut master_vec = Vec::new();
    let mut line_count = 1;
    let mut header = "".to_string();
    let no_hits_num: usize = colors_accession.len();
    let batch = b * 4;
    let num_files = filenames.len();
    let f = File::open(&filenames[0]).expect("file not found");
    let d1 = MultiGzDecoder::new(f);
    let iter1 = io::BufReader::new(d1).lines();
    ThreadPoolBuilder::new()
        .num_threads(t)
        .build_global()
        .unwrap();
    let false_positive_p = false_prob_map(colors_accession, ref_kmers_in, bloom_size, num_hash);

    let my_bigsi: Arc<&std::collections::HashMap<usize, Vec<u8>>> = Arc::new(bigsi_map);
    let false_positive_p_arc: Arc<std::collections::HashMap<usize, f64>> =
        Arc::new(false_positive_p);
    for line in iter1 {
        let l = line.unwrap();
        if line_count % 4 == 1 {
            header = l.to_string();
        } else if line_count % 4 == 2 {
            vec.push(vec![header.to_owned(), l.to_owned()]);
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
                        (r[0].to_owned(), "too_short")
                    } else {
                        let map = if m == 0 {
                            kmer::kmerize_vector(vec![r[1].to_string()], k, d)
                        } else {
                            kmer::minimerize_vector(vec![r[1].to_string()], k, m, d)
                        };
                        let report =
                            search_index(&child_bigsi, &map, bloom_size, num_hash, no_hits_num);
                        if report.is_empty() {
                            (r[0].to_owned(), "no_hits")
                        } else {
                            (
                                r[0].to_owned(),
                                kmer_pol(
                                    &report,
                                    &map,
                                    &colors_accession,
                                    &child_fp,
                                    no_hits_num,
                                    fp_correct,
                                ),
                            )
                        }
                    }
                })
                .collect();
            master_vec.append(&mut c);
            eprint!("{} reads classified\r", master_vec.len());
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
                (r[0].to_owned(), "too_short")
            } else {
                let map = if m == 0 {
                    kmer::kmerize_vector(vec![r[1].to_string()], k, d)
                } else {
                    kmer::minimerize_vector(vec![r[1].to_string()], k, m, d)
                };
                let report = search_index(&child_bigsi, &map, bloom_size, num_hash, no_hits_num);
                if report.is_empty() {
                    (r[0].to_owned(), "no_hits")
                } else {
                    (
                        r[0].to_owned(),
                        kmer_pol(
                            &report,
                            &map,
                            &colors_accession,
                            &child_fp,
                            no_hits_num,
                            fp_correct,
                        ),
                    )
                }
            }
        })
        .collect();
    master_vec.append(&mut c);
    eprintln!("{} read pairs classified\n", master_vec.len());
    match search_time.elapsed() {
        Ok(elapsed) => {
            eprintln!(
                "Classified {} reads in {} seconds",
                master_vec.len(),
                elapsed.as_secs()
            );
        }
        Err(e) => {
            // an error occurred!
            eprintln!("Error: {:?}", e);
        }
    }
    let mut file =
        File::create(format!("{}_reads.txt", prefix)).expect("could not create outfile!");
    for id in master_vec {
        file.write_all(format!("{}\t {}\n", id.0, id.1).as_bytes())
            .expect("could not write results!");
        let count = tax_map.entry(id.1.to_string()).or_insert(0);
        *count += 1;
    }
    let mut count_file =
        File::create(format!("{}_counts.txt", prefix)).expect("could not create outfile!");
    for (k, v) in &tax_map {
        count_file
            .write_all(format!("{}\t {}\n", k, v).as_bytes())
            .expect("could not write count results!");
    }
    tax_map
}
