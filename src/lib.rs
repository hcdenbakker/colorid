extern crate flate2;
extern crate simple_bloom;
extern crate kmer_fa;
extern crate bit_vec;
extern crate murmurhash64;
extern crate bincode;
extern crate serde;
extern crate serde_json;
#[macro_use]
extern crate serde_derive;

use std::collections::HashMap;
use std::io;
use std::io::prelude::*;
use std::fs::File;
use simple_bloom::BloomFilter;
use murmurhash64::murmur_hash64a;
use bit_vec::BitVec;
use serde::Serialize;
use serde_json::Serializer;
use std::io::Write;
use std::error::Error;
use flate2::read::GzDecoder;



pub fn tab_to_map(filename: String) -> std::collections::HashMap<std::string::String, String> {
    let mut map = HashMap::new();
    let f = File::open(filename).expect("file not found");
    for line in io::BufReader::new(f).lines() {
        let l =  line.unwrap();
        let v: Vec<&str> = l.split("\t").collect();
        //println!("{} {}",v[0], v[1]);
        map.insert(String::from(v[0]), String::from(v[1]));
    }
    map
}

pub fn build_bigsi2(map: std::collections::HashMap<std::string::String, String>,bloom_size: usize,
     num_hash: usize, k_size: usize) -> (std::collections::HashMap<usize, bit_vec::BitVec>,
     std::collections::HashMap<usize, String >, std::collections::HashMap<String, usize >){
    //get a hash map with taxa and bit vector from bloom filter
    let mut bit_map = HashMap::new();
    let mut ref_kmer = HashMap::new();
    let mut accessions = Vec::new();
    for (accession, v) in &map{
        println!("Adding {} to BIGSI", accession);
        let vec = kmer_fa::read_fasta(v.to_string());
        let kmers = kmer_fa::kmerize_vector(vec, k_size);
        ref_kmer.insert(accession.to_string(), kmers.len());
        let mut filter = BloomFilter::new(bloom_size as usize, num_hash as usize);
        for (kmer, _ ) in &kmers {
            filter.insert(&kmer);
            }
        bit_map.insert(accession, filter.bits);
        }
    for (accession, _) in &map{
        accessions.push(accession);
        }
    //create hash table with colors for accessions
    let mut accession_colors = HashMap::new();
    let mut colors_accession = HashMap::new();
    for (c,s) in accessions.iter().enumerate(){
          accession_colors.insert(s.to_string(), c);
          colors_accession.insert(c, s.to_string());
    }
    let num_taxa = accessions.len();
    let mut bigsi_map = HashMap::new();
        for i in 0..bloom_size{
            let mut bitvec = bit_vec::BitVec::from_elem(num_taxa,false);
            for (t, s) in &bit_map{
                if s[i] == true{
                    bitvec.set(*accession_colors.get(&t.to_string()).unwrap(), true);
                }
            }
            bigsi_map.insert(i, bitvec);
        }       
   (bigsi_map, colors_accession, ref_kmer)
}


//there seems to be a ~18 s overhead to this function for larger BIGSIs
//add a hashmap with vectors containing k-mer coverages of uniquely placed k-mers to estimate
//coverage per taxon
pub fn search_bigsi(kmer_query: std::collections::HashMap<std::string::String, i32>, 
    bigsi_map: std::collections::HashMap<usize, bit_vec::BitVec>,
    colors_accession: std::collections::HashMap<usize, String >,
    bloom_size: usize, num_hash: usize, k_size: usize) -> (std::collections::HashMap<String, u64 >,
    std::collections::HashMap<String, Vec<f64> >) //hashmap with name and vector containing freqs taxon-specific kmers 
    {
    println!("Search! Collecting slices");
    let mut report = HashMap::new();
    let mut uniq_freqs = HashMap::new();
    for (k, _) in &kmer_query{
        let mut kmer_slices = Vec::new();
        for i in 0..num_hash{
            let bit_index  = murmur_hash64a(k.as_bytes(), i as u64) % bloom_size as u64;
            let bi = bit_index as usize;
            //if bigsy contains bit_index, safe bit_vec on that position to temp vec, else break
            if bigsi_map.contains_key(&bi){
                kmer_slices.push(bigsi_map.get(&bi).unwrap());
            }else{
                let count = report.entry(String::from("No hits!")).or_insert(0);
                *count += 1;
                break;
            }
        }
    //we have to deal with no hits!
        if kmer_slices.len() == num_hash as usize{
            let original_first = kmer_slices[0];
            let mut first = bit_vec::BitVec::from_elem(original_first.len(),false);
            for i in 0..first.len(){
                if original_first[i] == true{
                    first.set(i, true);
                }
            }
        for i in 1..num_hash{
            let j = i as usize;
            first.intersect(&kmer_slices[j]);
         }
        let mut hits = Vec::new();
        for i in 0..first.len(){
            if first[i] == true{
                hits.push(colors_accession.get(&i).unwrap());
            }
        }
        if hits.len() > 0{
            for h in &hits{
                let count = report.entry(h.to_string()).or_insert(0);
                *count += 1;
            }
            if hits.len() == 1{
                let key = hits[0];
                let value = *kmer_query.get(&k.to_string()).unwrap() as f64;
                uniq_freqs.entry(key.to_string()).or_insert(Vec::new()).push(value);

            }
        }else{
        let count = report.entry(String::from("No hits!")).or_insert(0);
        *count += 1;
        }
     }else{
            let count = report.entry(String::from("No hits!")).or_insert(0);
            *count += 1;
        }
   }
   (report, uniq_freqs)
}


/*for now naive per read IDer, for now just takes the taxon with the highest frequency
pub fn per_read_id(filename: String,
                   bigsi_map: std::collections::HashMap<usize, bit_vec::BitVec>,
                   colors_accession: std::collections::HashMap<usize, String >,
                   bloom_size: usize,
                   num_hash: usize,
                   k_size: usize) -> std::collections::HashMap<std::string::String, u32>{
   let mut f = File::open(filename).expect("file not found");
   let mut map = HashMap::new();
   let mut line_count = 1;
   let d = GzDecoder::new(f);
   for line in io::BufReader::new(d).lines() {
       let l =  line.unwrap();
       if line_count%4 == 2 {
           let bigsi = bigsi_map.clone();
           let colors = colors_accession.clone();
           let read_vec = vec![l];
           let read_kmers = kmer_fa::kmerize_vector(read_vec, k_size);
           //we cannot use this function...overhead to high, I assume because we have to load the
           //bigsi for every single iteration
           let read_report = search_bigsi(read_kmers, bigsi, colors, bloom_size, num_hash, k_size);
           let mut max = 0;
           let mut max_label = String::from("default!");
           for (k,v) in read_report{
               if v > max{
                   max = v;
                   max_label = k.to_owned();
               }
           }
           let count = map.entry(max_label.to_string()).or_insert(0);
           *count += 1;
           
       }
        line_count += 1;
        println!("{}", line_count);
    }

    map
}
*/
#[derive(Debug, Serialize, Deserialize)]
pub struct BigsyMap {
    pub bloom_size: usize,
    pub num_hash: usize,
    pub k_size: usize,
    pub colors: HashMap<usize, String >,
    pub map: HashMap<usize, Vec<u64>>,
    pub n_ref_kmers: HashMap<String, usize >
}

pub fn save_bigsi(bigsi: std::collections::HashMap<usize, bit_vec::BitVec>,
                  colors_accession: std::collections::HashMap<usize, String>,
                  n_ref_kmers_in: std::collections::HashMap<String, usize >,
                  bloom_size_in: usize, num_hash_in: usize, k_size_in: usize, path: &str)
    /*-> Result<(), Box<Error>>*/ {
    let mut bigsi_map = HashMap::new();
    for (k,v)  in bigsi{
        //create vector of integers from bitvec (v)
        let mut positions = Vec::new();
        for i in 0..v.len(){
            let j = i as u64;
            if v[i] == true{
            positions.push(j);
            }
        if positions.len() > 0 {
            bigsi_map.insert(k, positions.to_owned());
            }
         }
    }
    let mappy = BigsyMap {map: bigsi_map, colors: colors_accession, n_ref_kmers: n_ref_kmers_in, bloom_size: bloom_size_in, num_hash: num_hash_in, k_size: k_size_in};
    let serialized = serde_json::to_string(&mappy).unwrap();
    let mut writer = File::create(path).unwrap();
    writer.write_all(serialized.as_bytes());
    }
    
pub fn read_bigsi(path: &str) -> (std::collections::HashMap<usize, bit_vec::BitVec>,
                                  std::collections::HashMap<usize, String>,
                                  std::collections::HashMap<String, usize>,
                                  usize,
                                  usize,
                                  usize){
    let mut file = File::open(path).unwrap();
    let mut contents = String::new();
    file.read_to_string(&mut contents);
    let deserialized: BigsyMap = serde_json::from_str(&contents).unwrap();
    let mut bigsi_map = HashMap::new();
    //create a full hashmap with bitvectors (first naive implementation, takes care of 'empty'
    //bitslices...may not be necessary...
    //for i in 0..deserialized.bloom_size{
    //    bigsi_map.insert(i, bit_vec::BitVec::from_elem(deserialized.colors.len(),false));
    //}
    for (key, vector) in deserialized.map{
        let mut bitvector = bit_vec::BitVec::from_elem(deserialized.colors.len(),false);
        for j in vector{
            let z = j as usize;
            bitvector.set( z , true);
        }
        bigsi_map.insert(key, bitvector);
    }    
    (bigsi_map, deserialized.colors, deserialized.n_ref_kmers, deserialized.bloom_size, deserialized.num_hash, deserialized.k_size)
}

//https://codereview.stackexchange.com/questions/173338/calculate-mean-median-and-mode-in-rust 
pub fn mode(numbers: &Vec<f64>) -> usize {
    let mut occurrences = HashMap::new();

    for value in numbers {
        *occurrences.entry(*value as usize).or_insert(0) += 1;
    }

    occurrences.into_iter()
        .max_by_key(|&(_, count)| count)
        .map(|(val, _)| val)
        .expect("Cannot compute the mode of zero numbers")
}
