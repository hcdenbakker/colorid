//simple bloom filter implementation with publicly accessible bit vectors
extern crate bit_vec;
extern crate murmurhash64;

use bit_vec::BitVec;
use murmurhash64::murmur_hash64a;

pub struct BloomFilter {
    pub bits: BitVec, // Bit vector
    pub num_hashes: usize, // # of hashes needed
}

impl BloomFilter {
    pub fn new(m: usize, eta: usize) -> BloomFilter {
        BloomFilter {
            bits: BitVec::from_elem(m,false),
            num_hashes: eta,
        }

    }
    pub fn insert(&mut self, value: &str) {
    // Generate a bit index for each of the hash functions needed
    for i in 0..self.num_hashes {
       let bit_index = (murmur_hash64a(value.as_bytes(), i as u64) % (self.bits.len() as u64)) as u64;
            self.bits.set(bit_index as usize, true);
        }
    }
    
    pub fn contains(&self, value: &str) -> bool {
        for i in 0..self.num_hashes {
            let bit_index = (murmur_hash64a(value.as_bytes(), i as u64) % (self.bits.len() as u64)) as u64;

            if self.bits[bit_index as usize] == false {
                return false;
            }
        }
        return true;
    }

}

#[cfg(test)]
mod tests {
    use super::{BloomFilter};


    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
    #[test]
    fn set_filter() {
        let new_filter = BloomFilter::new(250000000, 4);
        println!("{}",new_filter.bits.len());
    }

    #[test]
    fn add_filter() {
        let mut new_filter = BloomFilter::new(250000000, 4);
        new_filter.insert("ATGC");
    }
    #[test]
    fn use_filter() {
        let mut new_filter = BloomFilter::new(250000000, 4);
        new_filter.insert("ATGC");
        assert_eq!(new_filter.contains("ATGT"), false);
        assert_eq!(new_filter.contains("ATGC"), true);
    }
}
