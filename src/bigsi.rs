use bincode::{deserialize_from, serialize};
use fnv;
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
    pub colors: fnv::FnvHashMap<usize, String>,
    pub map: fnv::FnvHashMap<usize, Vec<u8>>,
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

pub fn save_bigsi(
    path: &str,
    mappy: &BigsyMap
) {
    let serialized: Vec<u8> = serialize(mappy).unwrap();
    let mut writer = File::create(path).unwrap();
    writer
        .write_all(&serialized)
        .expect("problems preparing serialized data for writing");
}

pub fn read_bigsi(
    path: &str,
) -> BigsyMap 
{
    let mut reader = BufReader::new(File::open(path).expect("Can't open index!"));
    let deserialized: BigsyMap =
        deserialize_from(&mut reader).expect("can't deserialize");
    deserialized
}

pub fn save_bigsi_mini(
    path: &str,
    mappy: &BigsyMapMini
) {
    let serialized: Vec<u8> = serialize(&mappy).unwrap();
    let mut writer = File::create(path).unwrap();
    writer
        .write_all(&serialized)
        .expect("problems preparing serialized data for writing");
}

pub fn read_bigsi_mini(
    path: &str,
) -> BigsyMapMini
{
    let mut reader = BufReader::new(File::open(path).expect("Can't open index!"));
    let deserialized: BigsyMapMini =
        deserialize_from(&mut reader).expect("cant deserialize");
    deserialized
}
