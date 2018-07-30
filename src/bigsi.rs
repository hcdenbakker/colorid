extern crate serde;

use bincode::{serialize, deserialize_from, Infinite};
use std;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::io::Write;

#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct BigsyMap {
    pub bloom_size: usize,
    pub num_hash: usize,
    pub k_size: usize,
    pub colors: HashMap<usize, String>,
    pub map: HashMap<usize, Vec<u8>>,
    pub n_ref_kmers: HashMap<String, usize>,
}

#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct BigsyMapMini {
    pub bloom_size: usize,
    pub num_hash: usize,
    pub k_size: usize,
    pub m_size: usize,
    pub colors: HashMap<usize, String>,
    pub map: HashMap<usize, Vec<u8>>,
    pub n_ref_kmers: HashMap<String, usize>,
}

pub fn save_bigsi(
    bigsi_map: std::collections::HashMap<usize, Vec<u8>>,
    colors_accession: std::collections::HashMap<usize, String>,
    n_ref_kmers_in: std::collections::HashMap<String, usize>,
    bloom_size_in: usize,
    num_hash_in: usize,
    k_size_in: usize,
    path: &str,
) {
    let mappy = BigsyMap {
        map: bigsi_map,
        colors: colors_accession,
        n_ref_kmers: n_ref_kmers_in,
        bloom_size: bloom_size_in,
        num_hash: num_hash_in,
        k_size: k_size_in,
    };
    let serialized: Vec<u8> = serialize(&mappy, Infinite).unwrap();
    let mut writer = File::create(path).unwrap();
    writer
        .write_all(&serialized)
        .expect("problems preparing serialized data for writing");
}

pub fn read_bigsi(
    path: &str,
) -> (
    std::collections::HashMap<usize, Vec<u8>>,
    std::collections::HashMap<usize, String>,
    std::collections::HashMap<String, usize>,
    usize,
    usize,
    usize,
) {
    let mut reader = BufReader::new(File::open(path).expect("Can't open index!"));
    //let mut buffer = Vec::new();
    //reader.read_to_end(&mut buffer).expect("Can't read content");
    //let deserialized: BigsyMap = deserialize(&buffer[..]).expect("cant deserialize");
    let deserialized: BigsyMap = deserialize_from(&mut reader, Infinite).expect("cant deserialize");
    (
        deserialized.map,
        deserialized.colors,
        deserialized.n_ref_kmers,
        deserialized.bloom_size,
        deserialized.num_hash,
        deserialized.k_size,
    )
}

pub fn save_bigsi_mini(
    bigsi_map: std::collections::HashMap<usize, Vec<u8>>,
    colors_accession: std::collections::HashMap<usize, String>,
    n_ref_kmers_in: std::collections::HashMap<String, usize>,
    bloom_size_in: usize,
    num_hash_in: usize,
    k_size_in: usize,
    m_size_in: usize,
    path: &str,
) {
    let mappy = BigsyMapMini {
        map: bigsi_map,
        colors: colors_accession,
        n_ref_kmers: n_ref_kmers_in,
        bloom_size: bloom_size_in,
        num_hash: num_hash_in,
        k_size: k_size_in,
        m_size: m_size_in,
    };
    let serialized: Vec<u8> = serialize(&mappy, Infinite).unwrap();
    let mut writer = File::create(path).unwrap();
    writer
        .write_all(&serialized)
        .expect("problems preparing serialized data for writing");
}

pub fn read_bigsi_mini(
    path: &str,
) -> (
    std::collections::HashMap<usize, Vec<u8>>,
    std::collections::HashMap<usize, String>,
    std::collections::HashMap<String, usize>,
    usize,
    usize,
    usize,
    usize,
) {
    let mut reader = BufReader::new(File::open(path).expect("Can't open index!"));
    //let mut buffer = Vec::new();
    //reader.read_to_end(&mut buffer).expect("Can't read content");
    //let deserialized: BigsyMapMini = deserialize(&buffer[..]).expect("cant deserialize");
    let deserialized: BigsyMapMini = deserialize_from(&mut reader, Infinite).expect("cant deserialize");
    (
        deserialized.map,
        deserialized.colors,
        deserialized.n_ref_kmers,
        deserialized.bloom_size,
        deserialized.num_hash,
        deserialized.k_size,
        deserialized.m_size,
    )
}
