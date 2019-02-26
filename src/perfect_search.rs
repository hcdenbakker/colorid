use bit_vec::BitVec;
use fasthash;
use fnv;
use kmer;

pub fn batch_search(
    files: Vec<&str>,
    bigsi_map: &fnv::FnvHashMap<usize, BitVec>,
    colors_accession: &fnv::FnvHashMap<usize, String>,
    _n_ref_kmers: &fnv::FnvHashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
    k_size: usize,
    _cov: f64,
) {
    for file in files {
        //only fasta formatted file!
        eprintln!("Counting k-mers, this may take a while!");
        let vec_query = kmer::read_fasta(file.to_owned());
        let kmers_query = kmer::kmerize_vector(vec_query, k_size, 1);
        eprintln!("{} kmers in query", kmers_query.len());
        if kmers_query.len() == 0 {
            eprintln!("Warning! no kmers in query; maybe your kmer length is larger than your query length?");
        } else {
            let mut kmer_slices = Vec::new();
            for k in kmers_query.keys() {
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
            }
            if kmer_slices.len() < (num_hash * kmers_query.len()) {
                eprintln!("No perfect hits!");
            } else {
                //bit-wise AND
                let mut first = kmer_slices[0].to_owned();
                for i in 1..(num_hash * kmers_query.len()) {
                    let j = i as usize;
                    first.intersect(&kmer_slices[j]);
                }
                let mut hits = Vec::new();
                for i in 0..first.len() {
                    if first[i] {
                        hits.push(&colors_accession[&i]);
                    }
                }
                eprintln!("{} hits", hits.len());
                for h in &hits {
                    println!("{}\t{}\t{}\t1.00", file, h, kmers_query.len());
                }
            }
        }
    }
}

pub fn batch_search_mf(
    files: Vec<&str>,
    bigsi_map: &fnv::FnvHashMap<usize, BitVec>,
    colors_accession: &fnv::FnvHashMap<usize, String>,
    _n_ref_kmers: &fnv::FnvHashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
    k_size: usize,
    _cov: f64,
) {
    for file in files {
        let (labels, sequences) = super::kmer::read_fasta_mf(file.to_owned());
        for (i, label) in labels.iter().enumerate() {
            //only fasta formatted file!
            //eprintln!("Counting k-mers, this may take a while!");
            let kmers = super::kmer::kmerize_string(sequences[i].to_owned(), k_size);
            //let kmers_query = kmer_fa::clean_map(unfiltered, filter);
            match kmers{
                None    => println!("Warning! no kmers in query '{}'; maybe your kmer length is larger than your query length?", label),
                Some(kmers_query) =>
                {eprintln!("{} kmers in query", kmers_query.len());
                let mut kmer_slices = Vec::new();
                for k in kmers_query.keys() {
                    for i in 0..num_hash {
                        let bit_index = fasthash::xx::hash64_with_seed(&k.as_bytes(), i as u64)
                            % bloom_size as u64;
                        let bi = bit_index as usize;
                        if !bigsi_map.contains_key(&bi) {
                            break;
                        } else {
                            kmer_slices.push(&bigsi_map[&bi]);
                        }
                    }
                }
                if kmer_slices.len() < (num_hash * kmers_query.len()) {
                    eprintln!("No perfect hits!");
                } else {
                    //bit-wise AND
                    let mut first = kmer_slices[0].to_owned();
                    for i in 1..(num_hash * kmers_query.len()) {
                        let j = i as usize;
                        first.intersect(&kmer_slices[j]);
                    }
                    let mut hits = Vec::new();
                    for i in 0..first.len() {
                        if first[i] {
                            hits.push(&colors_accession[&i]);
                        }
                    }
                    eprintln!("{} hits", hits.len());
                    for h in &hits {
                        println!("{}\t{}\t{}\t1.00", label.to_string(), h, kmers_query.len());
                    }
                }
                },
            }
        }
    }
}
