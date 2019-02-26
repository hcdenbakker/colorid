use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::prelude::*;

pub fn tab_to_map(
    filename: String,
    query: &str,
    accept: &str,
) -> std::collections::HashMap<std::string::String, String> {
    let mut map = HashMap::new();
    let f = File::open(filename).expect("classification file not found");
    for line in io::BufReader::new(f).lines() {
        let l = line.unwrap();
        let v: Vec<&str> = l.split('\t').collect();
        let h: Vec<&str> = v[0].split(' ').collect();
        if v[1].contains(query)
        /* && (v[4] == accept)*/
        {
            map.insert(String::from(h[0]), String::from(v[1]));
        }
    }
    map
}

#[allow(unused_assignments)]
pub fn read_filter_pe(
    class_map: std::collections::HashMap<std::string::String, String>,
    filenames: Vec<&str>,
    query: &str,
    prefix: &str,
    exclude: bool,
) {
    let query_cleaned = query.replace(" ", "_");
    let mut line_count = 1;
    let f = File::open(&filenames[0]).expect("file not found");
    let f2 = File::open(&filenames[1]).expect("file not found");
    let d1 = MultiGzDecoder::new(f);
    let d2 = MultiGzDecoder::new(f2);
    let iter1 = io::BufReader::new(d1).lines();
    let mut iter2 = io::BufReader::new(d2).lines();
    let mut header1 = "".to_string();
    let mut header2 = "".to_string();
    let mut seq1 = "".to_string();
    let mut seq2 = "".to_string();
    let mut qual1 = "".to_string();
    let mut qual2 = "".to_string();
    let mut excluded = 0;
    let mut included = 0;
    let fq1 = File::create(format!("{}_{}_R1.fq.gz", prefix, query_cleaned))
        .expect("could not create R1!");
    let mut gz1 = GzEncoder::new(fq1, Compression::default());
    let fq2 = File::create(format!("{}_{}_R2.fq.gz", prefix, query_cleaned))
        .expect("could not create R2!");
    let mut gz2 = GzEncoder::new(fq2, Compression::default());
    for line in iter1 {
        let l = line.unwrap();
        let line2 = iter2.next();
        //maybe check here if headers come from same pair
        if line_count % 4 == 1 {
            match line2 {
                Some(h2) => {
                    header1 = l.to_string();
                    header2 = h2.unwrap().to_string();
                }
                None => break,
            };
        } else if line_count % 4 == 2 {
            match line2 {
                Some(l2) => {
                    seq1 = l.to_owned();
                    seq2 = l2.unwrap().to_owned();
                }
                None => break,
            };
        } else if line_count % 4 == 0 {
            match line2 {
                Some(l2) => {
                    qual1 = l.to_owned();
                    qual2 = l2.unwrap().to_owned();
                    let v: Vec<&str> = header1.split(' ').collect();
                    if exclude == true {
                        if !class_map.contains_key(v[0]) {
                            gz1.write_all(
                                format!("{}\n{}\n+\n{}\n", header1, seq1, qual1).as_bytes(),
                            )
                            .expect("could not write R1!");
                            gz2.write_all(
                                format!("{}\n{}\n+\n{}\n", header2, seq2, qual2).as_bytes(),
                            )
                            .expect("could not write R2!");
                            excluded += 1;
                        }
                    } else {
                        if class_map.contains_key(v[0]) {
                            gz1.write_all(
                                format!("{}\n{}\n+\n{}\n", header1, seq1, qual1).as_bytes(),
                            )
                            .expect("could not write R1!");
                            gz2.write_all(
                                format!("{}\n{}\n+\n{}\n", header2, seq2, qual2).as_bytes(),
                            )
                            .expect("could not write R2!");
                            included += 1;
                        }
                    }
                }
                None => break,
            };
        }
        line_count += 1;
    }
    gz1.finish().expect("Could not close new R1 file");
    gz2.finish().expect("Could not close new R2 file");
    if exclude == true {
        eprintln!(
            "Excluded {} read pairs  with classification containing '{}' from output files",
            excluded, query
        );
    } else {
        eprintln!(
            "Wrote {} read-pairs with classification containing '{}' to output files",
            included, query
        );
    }
}

//quicker to adjust the function to se
#[allow(unused_assignments)]
pub fn read_filter_se(
    class_map: std::collections::HashMap<std::string::String, String>,
    filenames: Vec<&str>,
    query: &str,
    prefix: &str,
    exclude: bool,
) {
    let query_cleaned = query.replace(" ", "_");
    let mut line_count = 1;
    let f = File::open(&filenames[0]).expect("file not found");
    let d1 = MultiGzDecoder::new(f);
    let iter1 = io::BufReader::new(d1).lines();
    let mut header1 = "".to_string();
    let mut seq1 = "".to_string();
    let mut qual1 = "".to_string();
    let mut excluded = 0;
    let mut included = 0;
    let fq1 =
        File::create(format!("{}_{}.fq.gz", prefix, query_cleaned)).expect("could not create R1!");
    let mut gz1 = GzEncoder::new(fq1, Compression::default());
    for line in iter1 {
        let l = line.unwrap();
        if line_count % 4 == 1 {
            header1 = l.to_string();
        } else if line_count % 4 == 2 {
            seq1 = l.to_owned();
        } else if line_count % 4 == 0 {
            qual1 = l.to_owned();
            let v: Vec<&str> = header1.split(' ').collect();
            if exclude == true {
                if !class_map.contains_key(v[0]) {
                    gz1.write_all(format!("{}\n{}\n+\n{}\n", header1, seq1, qual1).as_bytes())
                        .expect("Could not write forward read(-s) to file");
                    excluded += 1;
                }
            } else {
                if class_map.contains_key(v[0]) {
                    gz1.write_all(format!("{}\n{}\n+\n{}\n", header1, seq1, qual1).as_bytes())
                        .expect("Could not write reverse read(-s) to file");
                    included += 1;
                }
            }
        }
        line_count += 1;
    }
    gz1.finish().expect("Could not close new read file");
    if exclude == true {
        eprintln!(
            "Excluded {} read pairs  with classification containing '{}' from output files",
            excluded, query
        );
    } else {
        eprintln!(
            "Wrote {} read-pairs with classification containing '{}' to output files",
            included, query
        );
    }
}
