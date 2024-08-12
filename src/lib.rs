extern crate bincode;
extern crate bit_vec;
extern crate xxh3;
extern crate flate2;
extern crate fnv;
extern crate probability;
extern crate rayon;
extern crate serde;
#[macro_use]
extern crate serde_derive;

pub mod build;

pub mod bigsi;

pub mod read_id_mt_pe;

pub mod read_id_batch;

pub mod perfect_search;

pub mod batch_search_pe;

pub mod kmer;

pub mod simple_bloom;

pub mod reports;

pub mod read_filter;

pub mod seq;
