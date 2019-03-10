use bigsi;
use build;
use read_id_mt_pe;
use reports;
use std::time::SystemTime;

pub fn read_id_batch(
    batch_samples: &str,
    index: &str,
    threads: usize,
    down_sample: usize,
    fp_correct: f64,
    batch: usize,
    quality: u8,
    bitvector_sample: usize,
    high_mem_load: bool,
    tag: &str,
) {
    //tab to map to read batch
    let batch_map = build::tab_to_map(batch_samples.to_string());
    let bigsi_time = SystemTime::now();
    if index.ends_with(".mxi") {
        //let metadata = fs::metadata(&index).expect("Can't read metadata index!");
        let bigsi = if high_mem_load {
            bigsi::read_bigsi_mini_highmem(index)
        } else {
            bigsi::read_bigsi_mini(index)
        };
        match bigsi_time.elapsed() {
            Ok(elapsed) => {
                eprintln!("Index loaded in {} seconds", elapsed.as_secs());
            }
            Err(e) => {
                // an error occurred!
                eprintln!("Error: {:?}", e);
            }
        }
        //start iterating through files here
        for (accession, fq_vec) in &batch_map {
            eprintln!("Classifying {}", accession);
            let prefix = format!("{}_{}", accession, tag);
            let mut fq: Vec<_> = vec![];
            fq = fq_vec.iter().map(|r|{&r[..]}).collect();
            if fq[0].ends_with(".gz") {
                if fq.len() > 1 {
                    read_id_mt_pe::per_read_stream_pe(
                        fq,
                        &bigsi.map,
                        &bigsi.colors,
                        &bigsi.n_ref_kmers,
                        bigsi.bloom_size,
                        bigsi.num_hash,
                        bigsi.k_size,
                        bigsi.m_size,
                        threads,
                        down_sample,
                        fp_correct,
                        batch,
                        &prefix,
                        quality,
                        bitvector_sample,
                    )
                } else {
                    read_id_mt_pe::per_read_stream_se(
                        fq,
                        &bigsi.map,
                        &bigsi.colors,
                        &bigsi.n_ref_kmers,
                        bigsi.bloom_size,
                        bigsi.num_hash,
                        bigsi.k_size,
                        bigsi.m_size,
                        threads,
                        down_sample,
                        fp_correct,
                        batch,
                        &prefix,
                        quality,
                        bitvector_sample,
                    )
                };
            } else {
                read_id_mt_pe::stream_fasta(
                    fq,
                    &bigsi.map,
                    &bigsi.colors,
                    &bigsi.n_ref_kmers,
                    bigsi.bloom_size,
                    bigsi.num_hash,
                    bigsi.k_size,
                    bigsi.m_size,
                    threads,
                    down_sample,
                    fp_correct,
                    batch,
                    &prefix,
                    bitvector_sample,
                );
            }
            reports::read_counts_five_fields(prefix.to_owned() + "_reads.txt", &prefix);
        }
    } else {
        let bigsi = if high_mem_load {
            bigsi::read_bigsi_highmem(index)
        } else {
            bigsi::read_bigsi(index)
        };
        match bigsi_time.elapsed() {
            Ok(elapsed) => {
                eprintln!("Index loaded in {} seconds", elapsed.as_secs());
            }
            Err(e) => {
                // an error occurred!
                eprintln!("Error: {:?}", e);
            }
        }
        for (accession, fq_vec) in &batch_map {
            eprintln!("Classifying {}", accession);
            let prefix = format!("{}_{}", accession, tag);
            let mut fq: Vec<_> = vec![];
            fq = fq_vec.iter().map(|r|{&r[..]}).collect();
            if fq[0].ends_with(".gz") {
                if fq.len() > 1 {
                    read_id_mt_pe::per_read_stream_pe(
                        fq,
                        &bigsi.map,
                        &bigsi.colors,
                        &bigsi.n_ref_kmers,
                        bigsi.bloom_size,
                        bigsi.num_hash,
                        bigsi.k_size,
                        0,
                        threads,
                        down_sample,
                        fp_correct,
                        batch,
                        &prefix,
                        quality,
                        bitvector_sample,
                    )
                } else {
                    read_id_mt_pe::per_read_stream_se(
                        fq,
                        &bigsi.map,
                        &bigsi.colors,
                        &bigsi.n_ref_kmers,
                        bigsi.bloom_size,
                        bigsi.num_hash,
                        bigsi.k_size,
                        0,
                        threads,
                        down_sample,
                        fp_correct,
                        batch,
                        &prefix,
                        quality,
                        bitvector_sample,
                    )
                };
            } else {
                read_id_mt_pe::stream_fasta(
                    fq,
                    &bigsi.map,
                    &bigsi.colors,
                    &bigsi.n_ref_kmers,
                    bigsi.bloom_size,
                    bigsi.num_hash,
                    bigsi.k_size,
                    0,
                    threads,
                    down_sample,
                    fp_correct,
                    batch,
                    &prefix,
                    bitvector_sample,
                );
            }
            reports::read_counts_five_fields(prefix.to_owned() + "_reads.txt", &prefix);
        }
    }
}
