extern crate bigs_id;
#[macro_use]
extern crate clap;
extern crate kmer_fa;

use clap::{App, AppSettings, Arg, SubCommand};
use std::time::SystemTime;

fn main() {
    let matches = App::new("bigs_id")
        .version("0.4")
        .author("Henk C. den Bakker <henkcdenbakker@gmail.com>")
        .about("BIGSI based taxonomic ID of sequence data")
        .setting(AppSettings::ArgRequiredElseHelp)
        .subcommand(
            SubCommand::with_name("build")
                .about("builds a bigsi")
                .version("0.1")
                .author("Henk C. den Bakker <henkcdenbakker@gmail.com>")
                .setting(AppSettings::ArgRequiredElseHelp)
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
                              -c, --compressed 'if set to 'true', will create a compressed index (default: false)'
                              -t, --threads 'sets number of threads, if set to 0, takes as many threads as it can get, default 1'")
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
                    Arg::with_name("threads")
                        .help("number of threads to use, if not set one thread will be used")
                        .required(false)
                        .short("t")
                        .takes_value(true)
                        .long("threads"),
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
                .about("does a bigsi search on one or more fasta/fastq.gz files")
                .version("0.1")
                .author("Henk C. den Bakker <henkcdenbakker@gmail.com>")
                .setting(AppSettings::ArgRequiredElseHelp)
                .arg(
                    Arg::with_name("bigsi")
                        .help("Sets the name of the index file for search")
                        .short("b")
                        .long("bigsi")
                        .required(true)
                        .takes_value(true),
                )
                .help(
                              "                              -b, --bigsi=[FILE] 'Sets the name of the index file for search'
                              -q, --query      'one or more fastq.gz or fasta formatted files to be queried'
                              -r, --reverse    'one or more fastq.gz (not fasta!) formatted files, supposed to be reverse reads of a paired end run'  
                              -f, --filter     'Sets threshold to filter k-mers by frequency'
                              -p, --p_shared        'minimum proportion of kmers shared with reference'
                              -c, --compressed 'If set to 'true', will assume compressed index (default: false)'
                              -g, If set('-g'), the proportion of kmers from the query matching the entries in the index will be reported'
                              -s, If set('-s'), the 'speedy' 'perfect match' alforithm will be performed'")
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
                    Arg::with_name("reverse")
                        .help("reverse file(-s)fastq.gz")
                        .required(false)
                        .min_values(1)
                        .default_value("none")
                        .short("r")
                        .takes_value(true)
                        .long("reverse"),
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
                        .help("If set('-g'), the proportion of kmers from the query matching the entries in the index will be reported")
                        .required(false)
                        .short("g")
                        .takes_value(false)
                        .long("gene_search"),
                )
                .arg(
                    Arg::with_name("perfect_search")
                        .help("If ('-s') is set, the fast 'perfect match' algorithm will be used")
                        .required(false)
                        .short("s")
                        .takes_value(false)
                        .long("perfect_search"),
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
                .setting(AppSettings::ArgRequiredElseHelp)
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
        .subcommand(
            SubCommand::with_name("read_id")
                .about("id's reads")
                .version("0.2")
                .author("Henk den Bakker. <henkcdenbakker@gmail.com>")
                .setting(AppSettings::ArgRequiredElseHelp)
                .arg(
                    Arg::with_name("bigsi")
                        .short("b")
                        .long("bigsi")
                        .required(true)
                        .takes_value(true)
                        .help("index for which info is requested"),
                        )
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
                    Arg::with_name("compressed")
                        .help("If set to 'true', it is assumed a gz compressed index is used")
                        .required(false)
                        .short("c")
                        .takes_value(true)
                        .long("compressed"),
                )
                .arg(
                    Arg::with_name("threads")
                        .help("number of threads to use, if not set the maximum available threads will be used")
                        .required(false)
                        .short("t")
                        .takes_value(true)
                        .long("threads"),
                )
                .arg(
                    Arg::with_name("down_sample")
                        .help("down-sample k-mers used for read classification, default 1; increases speed at cost of decreased sensitivity ")
                        .required(false)
                        .short("d")
                        .takes_value(true)
                        .long("down_sample"),
                )
                .arg(
                    Arg::with_name("fp_correct")
                        .help("parameter to correct for false positives, default 0.001, maybe increased for larger searches")
                        .required(false)
                        .short("p")
                        .takes_value(true)
                        .long("fp_correct"),
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
        let kmer = value_t!(matches, "k-mer_size", usize).unwrap_or(31);
        let bloom = value_t!(matches, "length_bloom", usize).unwrap_or(50_000_000);
        let hashes = value_t!(matches, "num_hashes", usize).unwrap_or(4);
        let threads = value_t!(matches, "threads", usize).unwrap_or(1);
        let compressed = value_t!(matches, "compressed", bool).unwrap_or(false);
        let map = bigs_id::tab_to_map(matches.value_of("ref_file").unwrap().to_string());
        let (bigsi_map, colors_accession, n_ref_kmers) = if threads == 1 {
            bigs_id::build_bigsi2(&map, bloom, hashes, kmer)
        } else {
            bigs_id::build_mt::build_bigsi(&map, bloom, hashes, kmer, threads)
        };
        println!("Saving BIGSI to file.");
        if !compressed {
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
        let files1: Vec<_> = matches.values_of("query").unwrap().collect();
        let files2 = if matches.value_of("reverse").unwrap() == "none"{
            vec![]
        }else{
            matches.values_of("reverse").unwrap().collect()
        };
        //let files2: Vec<_> = matches.values_of("reverse").unwrap().collect();
        let filter = value_t!(matches, "filter", isize).unwrap_or(-1);
        let cov = value_t!(matches, "shared_kmers", f64).unwrap_or(0.35);
        let gene_search = matches.is_present("gene_search");
        let perfect_search = matches.is_present("perfect_search");
        let compressed = value_t!(matches, "compressed", bool).unwrap_or(false);
        let bigsi_time = SystemTime::now();
        eprintln!("Loading index");
        let (bigsi_map, colors_accession, n_ref_kmers, bloom_size, num_hash, k_size) =
            if !compressed {
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
        if perfect_search {
            bigs_id::perfect_search::batch_search(
                files1,
                &bigsi_map,
                &colors_accession,
                &n_ref_kmers,
                bloom_size,
                num_hash,
                k_size,
                cov,
            )
        } else {
            bigs_id::batch_search_pe::batch_search(
                files1,
                files2,
                &bigsi_map,
                &colors_accession,
                &n_ref_kmers,
                bloom_size,
                num_hash,
                k_size,
                filter,
                cov,
                gene_search,
            )
        }
    }

    if let Some(matches) = matches.subcommand_matches("info") {
        let compressed = value_t!(matches, "compressed", bool).unwrap_or(false);
        let bigsi_time = SystemTime::now();
        eprintln!("Loading index");
        let (_bigsi_map, colors_accession, _n_ref_kmers, bloom_size, num_hash, k_size) =
            if !compressed {
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
        println!(
            "BIGSI parameters:\nBloomfilter-size: {}\nNumber of hashes: {}\nK-mer size: {}",
            bloom_size, num_hash, k_size
        );
        println!("Number of accessions in index: {}", colors_accession.len());
        let mut accessions = Vec::new();
        for (_k, v) in colors_accession {
            accessions.push(v);
        }
        accessions.sort_by(|a, b| a.cmp(b));
        for a in accessions {
            println!("{}", a);
        }
        //let accessions_sorted = accessions.sort_by(|a, b| b.cmp(a));
        //println!("{}", accessions_sorted.len());
    }
    if let Some(matches) = matches.subcommand_matches("read_id") {
        let bigsi_time = SystemTime::now();
        //let fq = matches.value_of("query").unwrap();
        let fq: Vec<_> = matches.values_of("query").unwrap().collect();
        let compressed = value_t!(matches, "compressed", bool).unwrap_or(false);
        let threads = value_t!(matches, "threads", usize).unwrap_or(0);
        let down_sample = value_t!(matches, "down_sample", usize).unwrap_or(1);
        let fp_correct = value_t!(matches, "fp_correct", f64).unwrap_or(0.001);
        eprintln!("Loading index");
        let (bigsi_map, colors_accession, n_ref_kmers, bloom_size, num_hash, k_size) =
            if !compressed {
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
        let tax_map = bigs_id::read_id_mt_pe::per_read_search(
            fq,
            bigsi_map,
            &colors_accession,
            &n_ref_kmers,
            bloom_size,
            num_hash,
            k_size,
            threads,
            down_sample,
            fp_correct,
        );
        for (k, v) in tax_map {
            println!("{}: {}", k, v);
        }
    }
}
