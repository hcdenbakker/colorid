pub struct Fasta {
    pub id: String,  //id including >
    pub seq: String, // sequence
}

impl Fasta {
    pub fn new() -> Fasta {
        Fasta {
            id: String::new(),
            seq: String::new(),
        }
    }
}

pub struct Fastq {
    pub id: String,    //id including >
    pub seq1: String,  // sequence 1
    pub qual1: String, // qual seq 1
    pub seq2: String,
    pub qual2: String,
}

impl Fastq {
    pub fn new() -> Fastq {
        Fastq {
            id: String::new(),
            seq1: String::new(),
            qual1: String::new(),
            seq2: String::new(),
            qual2: String::new(),
        }
    }
}

//from L. Katz fasten! make a method for Fastq
pub fn qual_mask(seq: &String, qual: &String, max_quality_offset: u8) -> String {
    if max_quality_offset == 0 {
        seq.to_string()
    } else {
        let max_quality: u8 = max_quality_offset + 33;
        let mut seq_chars = seq.chars();
        let mut new_seq = String::new();
        for q in qual.chars() {
            let current_nt = seq_chars
                .next()
                .expect("ERROR: could not get the next nt in the sequence");
            let phred = q as u8;
            if phred < max_quality {
                new_seq.push('N');
            } else {
                new_seq.push(current_nt);
            }
        }
        new_seq
    }
}

//from needletail https://github.com/onecodex/needletail/blob/master/src/kmer.rs
pub fn is_good_base(chr: u8) -> bool {
    match chr as char {
        'a' | 'c' | 'g' | 't' | 'A' | 'C' | 'G' | 'T' => true,
        _ => false,
    }
}

pub fn has_no_n<'a>(seq: &'a [u8]) -> bool {
    //! Determines if a sequence has any non-primary four bases
    //! characters in it
    seq.iter().all(|n| is_good_base(*n))
}

#[test]
fn can_detect_no_n() {
    assert!(has_no_n(b"AAGT"));
    assert!(!has_no_n(b"NAGT"));
}
