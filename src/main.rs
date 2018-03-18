extern crate bigs_id;
#[macro_use]
extern crate clap;
extern crate kmer_fa;

use clap::{App, Arg, SubCommand};
use std::time::SystemTime;

fn main() {
    let matches = App::new("bigsID")
        .version("0.2")
        .author("Henk C. den Bakker <henkcdenbakker@gmail.com>")
        .about("BIGSI based taxonomic ID of sequence data")
        .subcommand(
            SubCommand::with_name("build")
                .about("builds a bigsi")
                .version("0.1")
                .author("Henk C. den Bakker <henkcdenbakker@gmail.com>")
                .arg(
                    Arg::with_name("bigsi")
                        .short("b")
                        .long("bigsi")
                        .required(true)
                        .takes_value(true),
                )
                .help(
                              "                              -b, --bigsi=[FILE] 'Sets the name of the index file'
                              -r, --refs      'two column tab delimited file, first column taxon name, second column file  with sequence data plus path'
                              -k, --kmer      'sets kmer size to use for index'
                              -n, --num_hashes  'number of hashes to use for bloom filters'
                              -s, --bloom     'size bloom filter'
                              -c, --compressed 'if set to 'true', will create a compressed index (default: false)'")
                .arg(
                    Arg::with_name("ref_file")
                        .help("Sets the input file to use")
                        .required(true)
                        .short("r")
                        .takes_value(true)
                        .long("refs"),
                )
                .arg(
                    Arg::with_name("k-mer_size")
                        .help("Sets k-mer size")
                        .required(true)
                        .short("k")
                        .takes_value(true)
                        .long("kmer"),
                )
                .arg(
                    Arg::with_name("num_hashes")
                        .help("Sets number of hashes for bloom filter")
                        .required(true)
                        .short("n")
                        .takes_value(true)
                        .long("num_hashes"),
                )
                .arg(
                    Arg::with_name("length_bloom")
                        .help("Sets the length of the bloom filter")
                        .required(true)
                        .short("s")
                        .takes_value(true)
                        .long("bloom"),
                )
                .arg(
                    Arg::with_name("compressed")
                        .help("If set to 'true', a gz compressed index is build")
                        .required(false)
                        .short("c")
                        .takes_value(true)
                        .long("compressed"),
                ),
        )
        .subcommand(
            SubCommand::with_name("search")
                .about("does a bigsi search")
                .version("0.1")
                .author("Henk C. den Bakker <henkcdenbakker@gmail.com>")
                .arg(
                    Arg::with_name("bigsi")
                        .short("b")
                        .long("bigsi")
                        .required(true)
                        .takes_value(true),
                )
                .help(
                              "                              -b, --bigsi=[FILE] 'Sets the name of the index file for search'
                              -q, --query      'query file in fasta or fastq.gz format'
                              -f, --filter     'Sets threshold to filter k-mers by frequency'
                              -p, --p_shared        'minimum proportion of kmers shared with reference'
                              -c, --compressed 'if set to 'true', will assume compressed index (default: false)'")
                .arg(
                    Arg::with_name("query")
                        .help("query file (fastq.gz or fasta")
                        .required(true)
                        .short("q")
                        .takes_value(true)
                        .long("query"),
                )
                .arg(
                    Arg::with_name("filter")
                        .help("set minimum k-mer frequency ")
                        .required(false)
                        .short("f")
                        .takes_value(true)
                        .long("filter"),
                )
                .arg(
                    Arg::with_name("shared_kmers")
                        .help("set minimum proportion of shared k-mers with a reference")
                        .required(false)
                        .short("p")
                        .takes_value(true)
                        .long("p_shared"),
                )
                .arg(
                    Arg::with_name("gene_search")
                        .help("If set to 'true', the proportion of kmers from the query matching the entries in the index will be reported")
                        .required(false)
                        .short("g")
                        .takes_value(true)
                        .long("gene_search"),
                )
                .arg(
                    Arg::with_name("compressed")
                        .help("If set to 'true', it is assumed a gz compressed index is used")
                        .required(false)
                        .short("c")
                        .takes_value(true)
                        .long("compressed"),
                ),
        )
        .subcommand(
            SubCommand::with_name("batch_search")
                .about("does a bigsi search on a bunch of fastq.gz files")
                .version("0.1")
                .author("Henk C. den Bakker <henkcdenbakker@gmail.com>")
                .arg(
                    Arg::with_name("bigsi")
                        .short("b")
                        .long("bigsi")
                        .required(true)
                        .takes_value(true),
                )
                .help(
                              "                              -b, --bigsi=[FILE] 'Sets the name of the index file for search'
                              -q, --query      'one or more fastq.gz formatted files to be queried'
                              -f, --filter     'Sets threshold to filter k-mers by frequency'
                              -p, --p_shared        'minimum proportion of kmers shared with reference'
                              -c, --compressed 'if set to 'true', will assume compressed index (default: false)'")
                .arg(
                    Arg::with_name("query")
                        .help("query file(-s)fastq.gz")
                        .required(true)
                        .min_values(1)
                        .short("q")
                        .takes_value(true)
                        .long("query"),
                )
                .arg(
                    Arg::with_name("filter")
                        .help("set minimum k-mer frequency ")
                        .required(false)
                        .short("f")
                        .takes_value(true)
                        .long("filter"),
                )
                .arg(
                    Arg::with_name("shared_kmers")
                        .help("set minimum proportion of shared k-mers with a reference")
                        .required(false)
                        .short("p")
                        .takes_value(true)
                        .long("p_shared"),
                )
                .arg(
                    Arg::with_name("gene_search")
                        .help("If set to 'true', the proportion of kmers from the query matching the entries in the index will be reported")
                        .required(false)
                        .short("g")
                        .takes_value(true)
                        .long("gene_search"),
                )
                .arg(
                    Arg::with_name("compressed")
                        .help("If set to 'true', it is assumed a gz compressed index is used")
                        .required(false)
                        .short("c")
                        .takes_value(true)
                        .long("compressed"),
                ),
        )
        /*.subcommand(
            SubCommand::with_name("test")
                .about("controls testing features")
                .version("0.1")
                .author("Someone E. <someone_else@other.com>")
                .arg(
                    Arg::with_name("debug")
                        .short("d")
                        .help("print debug information verbosely"),
                ),
        )*/
        .subcommand(
            SubCommand::with_name("info")
                .about("dumps index parameters and accessions")
                .version("0.1")
                .author("Henk den Bakker. <henkcdenbakker@gmail.com>")
                .arg(
                    Arg::with_name("bigsi")
                        .short("b")
                        .long("bigsi")
                        .required(true)
                        .takes_value(true)
                        .help("index for which info is requested"),
                        )
                .arg(
                    Arg::with_name("compressed")
                        .help("If set to 'true', it is assumed a gz compressed index is used")
                        .required(false)
                        .short("c")
                        .takes_value(true)
                        .long("compressed"),
                ),
        )
        .get_matches();

    if let Some(matches) = matches.subcommand_matches("build") {
        println!(" Ref_file : {}", matches.value_of("ref_file").unwrap());
        println!(" Bigsi file : {}", matches.value_of("bigsi").unwrap());
        println!("K-mer size: {}", matches.value_of("k-mer_size").unwrap());
        println!(
            "Bloom filter parameters: num hashes {}, filter size {}",
            matches.value_of("num_hashes").unwrap(),
            matches.value_of("length_bloom").unwrap()
        );
        println!("Building BIGSI.");
        let kmer = value_t!(matches, "k-mer_size", usize).unwrap_or(31);
        let bloom = value_t!(matches, "length_bloom", usize).unwrap_or(50000000);
        let hashes = value_t!(matches, "num_hashes", usize).unwrap_or(4);
        let compressed = value_t!(matches, "compressed", bool).unwrap_or(false);
        let map = bigs_id::tab_to_map(matches.value_of("ref_file").unwrap().to_string());
        let (bigsi_map, colors_accession, n_ref_kmers) =
            bigs_id::build_bigsi2(map, bloom, hashes, kmer);
        println!("Saving BIGSI to file.");
        if compressed == false {
            bigs_id::save_bigsi(
                bigsi_map.to_owned(),
                colors_accession.to_owned(),
                n_ref_kmers.to_owned(),
                bloom,
                hashes,
                kmer,
                matches.value_of("bigsi").unwrap(),
            )
        } else {
            bigs_id::save_bigsi_gz(
                bigsi_map.to_owned(),
                colors_accession.to_owned(),
                n_ref_kmers.to_owned(),
                bloom,
                hashes,
                kmer,
                matches.value_of("bigsi").unwrap(),
            )
        };
    }
    if let Some(matches) = matches.subcommand_matches("search") {
        //read BIGSI
        let filter = value_t!(matches, "filter", i32).unwrap_or(0);
        let cov = value_t!(matches, "shared_kmers", f64).unwrap_or(0.35);
        let gene_search = value_t!(matches, "gene_search", bool).unwrap_or(false);
        let compressed = value_t!(matches, "compressed", bool).unwrap_or(false);
        let bigsi_time = SystemTime::now();
        eprintln!("Reading BIGSI");
        let (bigsi_map, colors_accession, n_ref_kmers, bloom_size, num_hash, k_size) =
            if compressed == false {
                bigs_id::read_bigsi(matches.value_of("bigsi").unwrap())
            } else {
                bigs_id::read_bigsi_gz(matches.value_of("bigsi").unwrap())
            };
        match bigsi_time.elapsed() {
            Ok(elapsed) => {
                eprintln!("Index read in {} seconds", elapsed.as_secs());
            }
            Err(e) => {
                // an error occurred!
                eprintln!("Error: {:?}", e);
            }
        }
        let quersy_in = matches.value_of("query").unwrap();
        if quersy_in.ends_with("gz") {
            let unfiltered = kmer_fa::kmers_from_fq(quersy_in.to_owned(), k_size);
            let kmers_query = kmer_fa::clean_map(unfiltered, filter);
            let num_kmers = kmers_query.len() as f64;
            println!("{} k-mers in query", num_kmers);
            let bigsi_search = SystemTime::now();
            let (report, freqs, _multi_freqs) = bigs_id::search_bigsi(
                kmers_query,
                bigsi_map,
                colors_accession,
                bloom_size,
                num_hash,
            );

            bigs_id::generate_report(report, freqs, n_ref_kmers, cov);
            match bigsi_search.elapsed() {
                Ok(elapsed) => {
                    eprintln!("Search completed in {} seconds.", elapsed.as_secs());
                }
                Err(e) => {
                    // an error occurred!
                    eprintln!("Error: {:?}", e);
                }
            }
        } else {
            let vec_query = kmer_fa::read_fasta(quersy_in.to_owned());
            let unfiltered = kmer_fa::kmerize_vector(vec_query, k_size);
            let kmers_query = kmer_fa::clean_map(unfiltered, filter);
            let num_kmers = kmers_query.len() as f64;
            println!("{} k-mers in query", num_kmers);
            let bigsi_search = SystemTime::now();
            let (report, freqs, _multi_freqs) = bigs_id::search_bigsi(
                kmers_query,
                bigsi_map,
                colors_accession,
                bloom_size,
                num_hash,
            );
            if gene_search == false {
                bigs_id::generate_report(report.to_owned(), freqs, n_ref_kmers, cov);
            } else {
                bigs_id::generate_report_gene(report, num_kmers as usize);
            }
            match bigsi_search.elapsed() {
                Ok(elapsed) => {
                    eprintln!("Search completed in {} seconds.", elapsed.as_secs());
                }
                Err(e) => {
                    // an error occurred!
                    eprintln!("Error: {:?}", e);
                }
            }
        }
    }
    if let Some(matches) = matches.subcommand_matches("batch_search") {
        let files: Vec<_> = matches.values_of("query").unwrap().collect();
        let filter = value_t!(matches, "filter", i32).unwrap_or(0);
        let cov = value_t!(matches, "shared_kmers", f64).unwrap_or(0.35);
        let gene_search = value_t!(matches, "gene_search", bool).unwrap_or(false);
        let compressed = value_t!(matches, "compressed", bool).unwrap_or(false);
        let bigsi_time = SystemTime::now();
        eprintln!("Loading index");
        let (bigsi_map, colors_accession, n_ref_kmers, bloom_size, num_hash, k_size) =
            if compressed == false {
                bigs_id::read_bigsi(matches.value_of("bigsi").unwrap())
            } else {
                bigs_id::read_bigsi_gz(matches.value_of("bigsi").unwrap())
            };
        match bigsi_time.elapsed() {
            Ok(elapsed) => {
                eprintln!("Index loaded in {} seconds", elapsed.as_secs());
            }
            Err(e) => {
                // an error occurred!
                eprintln!("Error: {:?}", e);
            }
        }
        bigs_id::batch_search(
            files,
            bigsi_map,
            colors_accession,
            n_ref_kmers,
            bloom_size,
            num_hash,
            k_size,
            filter,
            cov,
            gene_search,
        )
    }
    if let Some(matches) = matches.subcommand_matches("info") {
        let compressed = value_t!(matches, "compressed", bool).unwrap_or(false);
                let bigsi_time = SystemTime::now();
        eprintln!("Loading index");
        let (_bigsi_map, colors_accession, _n_ref_kmers, bloom_size, num_hash, k_size) =
            if compressed == false {
                bigs_id::read_bigsi(matches.value_of("bigsi").unwrap())
            } else {
                bigs_id::read_bigsi_gz(matches.value_of("bigsi").unwrap())
            };
        match bigsi_time.elapsed() {
            Ok(elapsed) => {
                eprintln!("Index loaded in {} seconds", elapsed.as_secs());
            }
            Err(e) => {
                // an error occurred!
                eprintln!("Error: {:?}", e);
            }
        }
        println!("BIGSI parameters:\nBloomfilter-size: {}\nNumber of hashes: {}\nK-mer size: {}", bloom_size, num_hash, k_size );
        println!("Number of accessions in index: {}", colors_accession.len());
        let mut accessions = Vec::new();
        for (_k, v) in colors_accession{
            accessions.push(v);
        }
        accessions.sort_by(|a, b| a.cmp(b));
        for a in accessions{
            println!("{}", a);
        }
        //let accessions_sorted = accessions.sort_by(|a, b| b.cmp(a));
        //println!("{}", accessions_sorted.len());
    }
}
