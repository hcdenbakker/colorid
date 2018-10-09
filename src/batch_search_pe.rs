use bit_vec::BitVec;
use fasthash;
use fnv;
use kmer;
use reports;
use std;
use std::collections::HashMap;
use std::time::SystemTime;

pub fn batch_search(
    files1: Vec<&str>,
    files2: Vec<&str>,
    bigsi_map: &fnv::FnvHashMap<usize, Vec<u8>>,
    colors_accession: &fnv::FnvHashMap<usize, String>,
    n_ref_kmers: &fnv::FnvHashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
    k_size: usize,
    filter: isize,
    cov: f64,
    gene_search: bool,
    qual_offset: u8,
) {
    for (i, file1) in files1.iter().enumerate() {
        if file1.ends_with("gz") {
            let unfiltered = if files2.len() == 0 {
                eprintln!("{}", file1);
                eprintln!("Counting k-mers, this may take a while!");
                kmer::kmers_from_fq_qual(file1.to_owned().to_string(), k_size, 1, qual_offset)
            } else {
                eprintln!("Paired end: {} {}", file1, files2[i]);
                eprintln!("Counting k-mers, this may take a while!");
                kmer::kmers_fq_pe_qual(vec![&file1, &files2[i]], k_size, 1, qual_offset)
            };
            let kmers_query = if filter < 0 {
                let cutoff = kmer::auto_cutoff(unfiltered.to_owned());
                kmer::clean_map(unfiltered, cutoff)
            } else {
                kmer::clean_map(unfiltered, filter as usize)
            };
            let num_kmers = kmers_query.len() as f64;
            eprintln!("{} k-mers in query", num_kmers);
            let bigsi_search = SystemTime::now();
            let mut report = HashMap::new();
            let mut uniq_freqs = HashMap::new();
            for k in kmers_query.keys() {
                let mut kmer_slices = Vec::new();
                for i in 0..num_hash {
                    let bit_index =
                        fasthash::xx::hash64_with_seed(&k.as_bytes(), i as u64) % bloom_size as u64;
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
                reports::generate_report(
                    file1,
                    report,
                    &uniq_freqs,
                    &n_ref_kmers,
                    num_kmers as usize,
                    cov,
                );
            } else {
                reports::generate_report_gene(file1, report, num_kmers as usize, cov);
            }
        } else {
            //else we assume it is a fasta formatted file!
            eprintln!("{}", file1);
            eprintln!("Counting k-mers, this may take a while!");
            let vec_query = kmer::read_fasta(file1.to_owned().to_string());
            let unfiltered = kmer::kmerize_vector(vec_query, k_size, 1);
            let kmers_query = if gene_search {
                kmer::clean_map(unfiltered, 0)
            } else if filter < 0 {
                eprintln!("no gene search");
                let cutoff = kmer::auto_cutoff(unfiltered.to_owned());
                kmer::clean_map(unfiltered, cutoff)
            } else {
                kmer::clean_map(unfiltered, filter as usize)
            };
            let num_kmers = kmers_query.len() as f64;
            eprintln!("{} k-mers in query", num_kmers);
            let mut report = HashMap::new();
            let mut uniq_freqs = HashMap::new();
            for k in kmers_query.keys() {
                let mut kmer_slices = Vec::new();
                for i in 0..num_hash {
                    let bit_index =
                        fasthash::xx::hash64_with_seed(&k.as_bytes(), i as u64) % bloom_size as u64;
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
                reports::generate_report(
                    file1,
                    report,
                    &uniq_freqs,
                    &n_ref_kmers,
                    num_kmers as usize,
                    cov,
                );
            } else {
                reports::generate_report_gene(file1, report, num_kmers as usize, cov);
            }
        }
    }
}
