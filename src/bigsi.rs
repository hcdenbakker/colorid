use bincode::{deserialize, deserialize_from, serialize};
use bit_vec::BitVec;
use fnv;
use std;
use std::collections::HashMap;
use std::fs;
use std::fs::File;
use std::fs::OpenOptions;
use std::io;
use std::io::prelude::*;
use std::io::BufReader;
use std::io::Write;
use std::path::Path;

#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct BigsyMap {
    pub bloom_size: usize,
    pub num_hash: usize,
    pub k_size: usize,
    pub colors: fnv::FnvHashMap<usize, String>,
    pub map: fnv::FnvHashMap<usize, Vec<u8>>,
    pub n_ref_kmers: fnv::FnvHashMap<String, usize>,
}

#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct BigsyMapNew {
    pub bloom_size: usize,
    pub num_hash: usize,
    pub k_size: usize,
    pub colors: fnv::FnvHashMap<usize, String>,
    pub map: fnv::FnvHashMap<usize, BitVec>,
    pub n_ref_kmers: fnv::FnvHashMap<String, usize>,
}

#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct BigsyMapMini {
    pub bloom_size: usize,
    pub num_hash: usize,
    pub k_size: usize,
    pub m_size: usize,
    pub colors: fnv::FnvHashMap<usize, String>,
    pub map: fnv::FnvHashMap<usize, Vec<u8>>,
    pub n_ref_kmers: fnv::FnvHashMap<String, usize>,
}

#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct BigsyMapMiniNew {
    pub bloom_size: usize,
    pub num_hash: usize,
    pub k_size: usize,
    pub m_size: usize,
    pub colors: fnv::FnvHashMap<usize, String>,
    pub map: fnv::FnvHashMap<usize, BitVec>,
    pub n_ref_kmers: fnv::FnvHashMap<String, usize>,
}

pub fn save_bigsi(path: &str, mappy: &BigsyMapNew) {
    let serialized: Vec<u8> = serialize(mappy).unwrap();
    let mut writer = File::create(path).unwrap();
    writer
        .write_all(&serialized)
        .expect("problems preparing serialized data for writing");
}

pub fn read_bigsi(path: &str) -> BigsyMapNew {
    let mut reader = BufReader::new(File::open(path).expect("Can't open index!"));
    let deserialized: BigsyMapNew = deserialize_from(&mut reader).expect("can't deserialize");
    deserialized
}

pub fn read_bigsi_highmem(path: &str) -> BigsyMapNew {
    let deserialized: BigsyMapNew =
        deserialize(&fs::read(path).expect("Can't open index!")).expect("cant deserialize");
    deserialized
}

pub fn save_bigsi_mini(path: &str, mappy: &BigsyMapMiniNew) {
    let serialized: Vec<u8> = serialize(&mappy).unwrap();
    let mut writer = File::create(path).unwrap();
    writer
        .write_all(&serialized)
        .expect("problems preparing serialized data for writing");
}

pub fn read_bigsi_mini(path: &str) -> BigsyMapMiniNew {
    let mut reader = BufReader::new(File::open(path).expect("Can't open index!"));
    let deserialized: BigsyMapMiniNew = deserialize_from(&mut reader).expect("cant deserialize");
    deserialized
}

pub fn read_bigsi_mini_highmem(path: &str) -> BigsyMapMiniNew {
    let deserialized: BigsyMapMiniNew =
        deserialize(&fs::read(path).expect("Can't open index!")).expect("cant deserialize");
    deserialized
}
