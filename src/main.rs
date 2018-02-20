#[macro_use]
extern crate clap;
extern crate bigs_id;
extern crate kmer_fa;

use clap::{Arg, App, SubCommand};
use std::time::{Duration, SystemTime};

fn main() {
    let matches = App::new("bigsID")
                          .version("1.0")
                          .author("Henk C. den Bakker <henkcdenbakker@gmail.com>")
                          .about("BIGSI based taxonomic ID of sequence data")
                          .subcommand(SubCommand::with_name("build")
                                      .about("builds a bigsi")
                                      .version("1.0")
                                      .author("Henk C. den Bakker <henkcdenbakker@gmail.com>")
                                      .arg(Arg::with_name("bigsi")
                                          .short("b")
                                          .long("bigsi")
                                          .required(true)
                                          .takes_value(true))
                                          .help("Sets file name to safe BIGSI")
                                      .arg(Arg::with_name("ref_file")
                                          .help("Sets the input file to use")
                                          .required(true)
                                          .short("r")
                                          .takes_value(true)
                                          .long("refs"))
                                      .arg(Arg::with_name("k-mer_size")
                                          .help("Sets k-mer size")
                                          .required(true)
                                          .short("k")
                                          .takes_value(true)
                                          .long("kmer"))
                                      .arg(Arg::with_name("num_hashes")
                                          .help("Sets number of hashes for bloom filter")
                                          .required(true)
                                          .short("n")
                                          .takes_value(true)
                                          .long("num_hashes"))
                                      .arg(Arg::with_name("length_bloom")
                                          .help("Sets the length of the bloom filter")
                                          .required(true)
                                          .short("s")
                                          .takes_value(true)
                                          .long("bloom"))
                                      )
                                .subcommand(SubCommand::with_name("search")
                                      .about("does a bigsi search")
                                      .version("1.0")
                                      .author("Henk C. den Bakker <henkcdenbakker@gmail.com>")
                                      .arg(Arg::with_name("bigsi")
                                          .short("b")
                                          .long("bigsi")
                                          .required(true)
                                          .takes_value(true))
                                          .help("Sets BIGSI index to use for search")
                                      .arg(Arg::with_name("query")
                                          .help("query file (fastq.gz or fasta")
                                          .required(true)
                                          .short("q")
                                          .takes_value(true)
                                          .long("query"))
                                      )            
                          .subcommand(SubCommand::with_name("test")
                                      .about("controls testing features")
                                      .version("1.3")
                                      .author("Someone E. <someone_else@other.com>")
                                      .arg(Arg::with_name("debug")
                                          .short("d")
                                          .help("print debug information verbosely")))
                          .get_matches();



    if let Some(matches) = matches.subcommand_matches("build") {
        println!(" Ref_file : {}", matches.value_of("ref_file").unwrap());
        println!(" Bigsi file : {}", matches.value_of("bigsi").unwrap());
        println!("K-mer size: {}", matches.value_of("k-mer_size").unwrap());
        println!("Bloom filter parameters: num hashes {}, filter size {}", matches.value_of("num_hashes").unwrap(), matches.value_of("length_bloom").unwrap());
        println!("Building BIGSI.");
        let kmer = value_t!(matches, "k-mer_size", usize).unwrap_or(31);
        let bloom = value_t!(matches, "length_bloom", usize).unwrap_or(50000000);
        let hashes = value_t!(matches, "num_hashes", usize).unwrap_or(4);
        let map = bigs_id::tab_to_map(matches.value_of("ref_file").unwrap().to_string());
        let (bigsi_map, colors_accession, n_ref_kmers) = bigs_id::build_bigsi2(map, bloom, hashes, kmer);
        println!("Saving BIGSI to file.");
        bigs_id::save_bigsi(bigsi_map.to_owned(), colors_accession.to_owned(), n_ref_kmers.to_owned(), bloom, hashes, kmer, matches.value_of("bigsi").unwrap());
   }
   if let Some(matches) = matches.subcommand_matches("search") {
       //read BIGSI
       let bigsi_time = SystemTime::now();
       println!("Reading BIGSI");
       let (bigsi_map, colors_accession, n_ref_kmers, bloom_size, num_hash, k_size) = bigs_id::read_bigsi(matches.value_of("bigsi").unwrap());
       match bigsi_time.elapsed() {
           Ok(elapsed) => {
               println!("{}", elapsed.as_secs());
       }
       Err(e) => {
           // an error occurred!
           println!("Error: {:?}", e);
       }
    }
    //now check if we can do a search
    let bigsi_search = SystemTime::now();
    let quersy_in = matches.value_of("query").unwrap();
    if quersy_in.ends_with("gz"){
        let unfiltered = kmer_fa::kmers_from_fq(quersy_in.to_owned(), k_size);
        let kmers_query = kmer_fa::clean_map(unfiltered, 0);
        let num_kmers = kmers_query.len() as f64;
        println!("{} k-mers in query", num_kmers);
        let bigsi_search = SystemTime::now();
        let (report, freqs) = bigs_id::search_bigsi(kmers_query, bigsi_map, colors_accession, bloom_size, num_hash, k_size);

        //let report = proto_bigsi::per_read_id(quersy_in.to_owned(), bigsi_map, colors_accession, bloom_size, num_hash, k_size);
        for (k,v) in report{
           let frequencies = freqs.get(&k.to_string());
           let mut mean: f64 = 0.0;
           let mut mode: usize = 0;
           match frequencies{
               Some(_x) => {
                   mean = frequencies.unwrap().iter().fold(0.0,|a, &b| a + b)/ frequencies.unwrap().len() as f64;
                   mode = bigs_id::mode(frequencies.unwrap());
               }
                   None => continue,
           }
           let n_kmers = n_ref_kmers.get(&k.to_string());
                match n_kmers {
                    // The division was valid
                    Some(_x) => {
                        let genome_cov = v as f64 / *n_kmers.unwrap() as f64;
                        if genome_cov > 0.35 {
                            println!("{}: {:.2} {:.2} {} {}", k, genome_cov, mean, mode,  frequencies.unwrap().len());
                        }
                    },
                    // The division was invalid
                    None    => continue,
                }
           //println!("{}: {:.2}", k, v as f64/num_kmers );
           //println!("{}: {}", k, v);
        }
       match bigsi_search.elapsed() {
       Ok(elapsed) => {
           println!("{}", elapsed.as_secs());
       }
       Err(e) => {
           // an error occurred!
           println!("Error: {:?}", e);
       }
    }
        }else{
        let vec_query = kmer_fa::read_fasta(quersy_in.to_owned());
        let kmers_query = kmer_fa::kmerize_vector(vec_query, k_size);
        let num_kmers = kmers_query.len() as f64;
        println!("{} k-mers in query", num_kmers);
        let bigsi_search = SystemTime::now();
        let (report, freqs) = bigs_id::search_bigsi(kmers_query, bigsi_map, colors_accession, bloom_size, num_hash, k_size);
            for (k,v) in report{
                let frequencies = freqs.get(&k.to_string());
                let mut mean: f64 = 0.0;
                match frequencies{
                    Some(_x) => {
                       mean = frequencies.unwrap().iter().fold(0.0,|a, &b| a + b)/ frequencies.unwrap().len() as f64;
                   }
                   None => continue,
                }
                let n_kmers = n_ref_kmers.get(&k.to_string());
                match n_kmers {
                    // The division was valid
                    Some(_x) => {
                        let genome_cov = v as f64 / *n_kmers.unwrap() as f64;
                        if genome_cov > 0.25 {
                            println!("{}: {:.2} {:.2} {}", k, genome_cov, mean, frequencies.unwrap().len());
                        }
                    },
                    // The division was invalid
                    None    => continue,
                }
                //println!("{}: {:.2}", k, v/ *n_ref_kmers.get(&k.to_string()).unwrap() as u64);
                //println!("{}: {} ", k, v);
            }
        match bigsi_search.elapsed() {
       Ok(elapsed) => {
           println!("{}", elapsed.as_secs());
       }
       Err(e) => {
           // an error occurred!
           println!("Error: {:?}", e);
       }
    }
    }   
   }
}
