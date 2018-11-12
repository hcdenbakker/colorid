extern crate colorid;
#[macro_use]
extern crate clap;

use clap::{App, AppSettings, Arg, SubCommand};
use colorid::bigsi;
use colorid::build;
use colorid::read_id_mt_pe::false_prob;
use std::alloc::System;
use std::fs;
use std::time::SystemTime;

#[global_allocator]
static GLOBAL: System = System;

fn main() {
    let matches = App::new("colorid")
        .version("0.1.1")
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
                              "                              -b, --bigsi=[FILE] 'Sets the prefix of the index file'
                              -r, --refs      'two column tab delimited file, first column taxon name, second column file  with sequence data plus path'
                              -k, --kmer      'sets kmer size to use for index'
                              -n, --num_hashes  'number of hashes to use for bloom filters'
                              -s, --bloom     'size bloom filter'
                              -m, --minimizer 'build index from minimizers'
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
                    Arg::with_name("minimizer")
                        .help("build index with minimizers, default length minimizer is 15, unless indicated otherwise")
                        .required(false)
                        .short("m")
                        .takes_value(false)
                        .long("minimizer"),
                )
                .arg(
                    Arg::with_name("value")
                        .help("sets length minimizer (default 15)")
                        .required(false)
                        .short("v")
                        .default_value("15")
                        .takes_value(true)
                        .long("value"),
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
                    Arg::with_name("quality")
                        .help("minimum phred score to keep basepairs within read (default 15)")
                        .required(false)
                        .short("Q")
                        .takes_value(true)
                        .long("quality"),
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
                /*.help(
                              "                              -b, --bigsi=[FILE] 'Sets the name of the index file for search'
                              -q, --query      'one or more fastq.gz or fasta formatted files to be queried'
                              -r, --reverse    'one or more fastq.gz (not fasta!) formatted files, supposed to be reverse reads of a paired end run'  
                              -f, --filter     'Sets threshold to filter k-mers by frequency'
                              -p, --p_shared        'minimum proportion of kmers shared with reference'
                              -g, If set('-g'), the proportion of kmers from the query matching the entries in the index will be reported'
                              -s, If set('-s'), the 'speedy' 'perfect match' alforithm will be performed'
                              -m, If set('-m'), each accession in a multifasta will betreated as a separate query, currently only with the -s option'")*/
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
                    Arg::with_name("multi_fasta")
                        .help("If set('-m'), each accession in a multifasta will betreated as a separate query, currently only with the -s option")
                        .required(false)
                        .short("m")
                        .takes_value(false)
                        .long("multi_fasta"),
                )
                .arg(
                    Arg::with_name("quality")
                        .help("minimum phred score to keep basepairs within read (default 15)")
                        .required(false)
                        .short("Q")
                        .takes_value(true)
                        .long("quality"),
                ),
        )
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
                        .help("index to be used for search"),
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
                    Arg::with_name("batch")
                        .help("Sets size of batch of reads to be processed in parallel, currently only implemented for minimizers")
                        .required(false)
                        .short("c")
                        .takes_value(true)
                        .long("batch"),
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
                    Arg::with_name("prefix")
                        .help("prefix for output file(-s)")
                        .required(true)
                        .short("n") //running out of options here!
                        .takes_value(true)
                        .long("prefix"),
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
                    Arg::with_name("high_mem_load")
                        .help("when this flag is set, a faster, but less memory efficient method to load the index is used ")
                        .required(false)
                        .short("H")
                        .takes_value(false)
                        .long("high_mem_load"),
                )
                .arg(
                    Arg::with_name("fp_correct")
                        .help("parameter to correct for false positives, default 3 (= 0.001), maybe increased for larger searches")
                        .required(false)
                        .short("p")
                        .takes_value(true)
                        .long("fp_correct"),
                )
                .arg(
                    Arg::with_name("quality")
                        .help("minimum phred score to keep basepairs within read (default 15)")
                        .required(false)
                        .short("Q")
                        .takes_value(true)
                        .long("quality"),
                        ),
        )
        .subcommand(
            SubCommand::with_name("read_filter")
                .about("filters reads")
                .version("0.1")
                .author("Henk den Bakker. <henkcdenbakker@gmail.com>")
                .setting(AppSettings::ArgRequiredElseHelp)
                .arg(
                    Arg::with_name("classification")
                        .short("c")
                        .long("classification")
                        .required(true)
                        .takes_value(true)
                        .help("tab delimited read classification file generated with the read_id subcommand"),
                        )
                .arg(
                    Arg::with_name("files")
                        .help("query file(-s)fastq.gz")
                        .required(true)
                        .min_values(1)
                        .short("f")
                        .takes_value(true)
                        .long("files"),
                )
                .arg(
                    Arg::with_name("taxon")
                        .help("taxon to be in- or excluded from the read file(-s)")
                        .required(true)
                        .short("t")
                        .takes_value(true)
                        .long("taxon"),
                        )
                .arg(
                    Arg::with_name("prefix")
                        .help("prefix for output file(-s)")
                        .required(true)
                        .short("p")
                        .takes_value(true)
                        .long("prefix"),
                        )
                .arg(
                    Arg::with_name("exclude")
                        .help("If set('-e or --exclude'), reads for which the classification contains the taxon name will be excluded")
                        .required(false)
                        .short("e")
                        .takes_value(false)
                        .long("exclude"),
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
        let quality = value_t!(matches, "quality", u8).unwrap_or(15);
        let minimizer = matches.is_present("minimizer");
        //hack to work around current clap bug with default values only being &str
        let minimizer_value: usize = matches.value_of("value").unwrap().parse::<usize>().unwrap();
        let map = build::tab_to_map(matches.value_of("ref_file").unwrap().to_string());
        if minimizer {
            println!("Build with minimizers, minimizer size: {}", minimizer_value);
            let (bigsi_map, colors_accession, n_ref_kmers) = if threads == 1 {
                build::build_single_mini(&map, bloom, hashes, kmer, minimizer_value)
            } else {
                build::build_multi_mini(&map, bloom, hashes, kmer, minimizer_value, threads)
            };
            println!("Saving BIGSI to file.");
            let index = bigsi::BigsyMapMiniNew {
                map: bigsi_map,
                colors: colors_accession,
                n_ref_kmers: n_ref_kmers,
                bloom_size: bloom,
                num_hash: hashes,
                k_size: kmer,
                m_size: minimizer_value,
            };
            bigsi::save_bigsi_mini(
                &(matches.value_of("bigsi").unwrap().to_owned() + ".mxi"),
                &index,
            );
        /*
            bigsi::save_bigsi_mini(
                bigsi_map.to_owned(),
                colors_accession.to_owned(),
                n_ref_kmers.to_owned(),
                bloom,
                hashes,
                kmer,
                minimizer_value,
                &(matches.value_of("bigsi").unwrap().to_owned() + ".mxi"),
            );*/
        } else {
            let (bigsi_map, colors_accession, n_ref_kmers) = if threads == 1 {
                build::build_single(&map, bloom, hashes, kmer)
            } else {
                build::build_multi(&map, bloom, hashes, kmer, threads)
            };
            let index = bigsi::BigsyMapNew {
                map: bigsi_map,
                colors: colors_accession,
                n_ref_kmers: n_ref_kmers,
                bloom_size: bloom,
                num_hash: hashes,
                k_size: kmer,
            };
            println!("Saving BIGSI to file.");
            bigsi::save_bigsi(
                &(matches.value_of("bigsi").unwrap().to_owned() + ".bxi"),
                &index,
            );
        }
    }
    if let Some(matches) = matches.subcommand_matches("search") {
        let files1: Vec<_> = matches.values_of("query").unwrap().collect();
        let files2 = if matches.value_of("reverse").unwrap() == "none" {
            vec![]
        } else {
            matches.values_of("reverse").unwrap().collect()
        };
        //let files2: Vec<_> = matches.values_of("reverse").unwrap().collect();
        let filter = value_t!(matches, "filter", isize).unwrap_or(-1);
        let cov = value_t!(matches, "shared_kmers", f64).unwrap_or(0.35);
        let gene_search = matches.is_present("gene_search");
        let perfect_search = matches.is_present("perfect_search");
        let multi_fasta = matches.is_present("multi_fasta");
        let quality = value_t!(matches, "quality", u8).unwrap_or(15);
        if matches.value_of("bigsi").unwrap().ends_with(".mxi") {
            eprintln!(
                "Error: An index with minimizers (.mxi) is used, but not available for this function"
            );
        } else {
            let bigsi_time = SystemTime::now();
            eprintln!("Loading index");
            let index = bigsi::read_bigsi(matches.value_of("bigsi").unwrap());
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
                if multi_fasta {
                    //make 'perfect' batch....
                    colorid::perfect_search::batch_search_mf(
                        files1,
                        &index.map,
                        &index.colors,
                        &index.n_ref_kmers,
                        index.bloom_size,
                        index.num_hash,
                        index.k_size,
                        cov,
                    )
                } else {
                    colorid::perfect_search::batch_search(
                        files1,
                        &index.map,
                        &index.colors,
                        &index.n_ref_kmers,
                        index.bloom_size,
                        index.num_hash,
                        index.k_size,
                        cov,
                    )
                }
            } else {
                colorid::batch_search_pe::batch_search(
                    files1,
                    files2,
                    &index.map,
                    &index.colors,
                    &index.n_ref_kmers,
                    index.bloom_size,
                    index.num_hash,
                    index.k_size,
                    filter,
                    cov,
                    gene_search,
                    quality,
                )
            }
        }
    }

    if let Some(matches) = matches.subcommand_matches("info") {
        let bigsi_time = SystemTime::now();
        eprintln!("Loading index");
        let index = matches.value_of("bigsi").unwrap();
        if index.ends_with(".mxi") {
            let bigsi = bigsi::read_bigsi_mini(index);
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
            "BIGSI parameters:\nBloomfilter-size: {}\nNumber of hashes: {}\nK-mer size: {}\n minimizer size: {}\n",
            bigsi.bloom_size, bigsi.num_hash, bigsi.k_size, bigsi.m_size
        );
            println!("Number of accessions in index: {}", bigsi.colors.len());
            let mut accessions = Vec::new();
            for (_k, v) in bigsi.colors {
                accessions.push(v);
            }
            accessions.sort_by(|a, b| a.cmp(b));
            for a in accessions {
                let k_size = bigsi.n_ref_kmers.get(&a).unwrap();
                println!(
                    "{} {} {:.3}",
                    a,
                    k_size,
                    false_prob(
                        bigsi.bloom_size as f64,
                        bigsi.num_hash as f64,
                        *k_size as f64
                    )
                );
            }
        } else {
            let bigsi = bigsi::read_bigsi(index);
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
                bigsi.bloom_size, bigsi.num_hash, bigsi.k_size
            );
            println!("Number of accessions in index: {}", bigsi.colors.len());
            let mut accessions = Vec::new();
            for (_k, v) in bigsi.colors {
                accessions.push(v);
            }
            accessions.sort_by(|a, b| a.cmp(b));
            for a in accessions {
                let k_size = bigsi.n_ref_kmers.get(&a).unwrap();
                println!(
                    "{} {} {:.3}",
                    a,
                    k_size,
                    false_prob(
                        bigsi.bloom_size as f64,
                        bigsi.num_hash as f64,
                        *k_size as f64
                    )
                );
            }
        }
    }
    if let Some(matches) = matches.subcommand_matches("read_id") {
        let bigsi_time = SystemTime::now();
        //let fq = matches.value_of("query").unwrap();
        let fq: Vec<_> = matches.values_of("query").unwrap().collect();
        let threads = value_t!(matches, "threads", usize).unwrap_or(0);
        let down_sample = value_t!(matches, "down_sample", usize).unwrap_or(1);
        let correct = value_t!(matches, "fp_correct", f64).unwrap_or(3.0);
        let fp_correct = 10f64.powf(-correct);
        let index = matches.value_of("bigsi").unwrap();
        let prefix = matches.value_of("prefix").unwrap();
        let quality = value_t!(matches, "quality", u8).unwrap_or(15);
        let batch = value_t!(matches, "batch", usize).unwrap_or(50000);
        let high_mem_load = matches.is_present("high_mem_load");

        if index.ends_with(".mxi") {
            //let metadata = fs::metadata(&index).expect("Can't read metadata index!");
            let bigsi = if high_mem_load {
                bigsi::read_bigsi_mini_highmem(index)
            } else {
                bigsi::read_bigsi_mini(index)
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
            if fq[0].ends_with(".gz") {
                if fq.len() > 1 {
                    colorid::read_id_mt_pe::per_read_stream_pe(
                        fq,
                        &bigsi.map,
                        &bigsi.colors,
                        &bigsi.n_ref_kmers,
                        bigsi.bloom_size,
                        bigsi.num_hash,
                        bigsi.k_size,
                        bigsi.m_size,
                        threads,
                        down_sample,
                        fp_correct,
                        batch,
                        prefix,
                        quality,
                    )
                } else {
                    colorid::read_id_mt_pe::per_read_stream_se(
                        fq,
                        &bigsi.map,
                        &bigsi.colors,
                        &bigsi.n_ref_kmers,
                        bigsi.bloom_size,
                        bigsi.num_hash,
                        bigsi.k_size,
                        bigsi.m_size,
                        threads,
                        down_sample,
                        fp_correct,
                        batch,
                        prefix,
                        quality,
                    )
                };
            } else {
                colorid::read_id_mt_pe::stream_fasta(
                    fq,
                    &bigsi.map,
                    &bigsi.colors,
                    &bigsi.n_ref_kmers,
                    bigsi.bloom_size,
                    bigsi.num_hash,
                    bigsi.k_size,
                    bigsi.m_size,
                    threads,
                    down_sample,
                    fp_correct,
                    batch,
                    prefix,
                );
            }
        } else {
            let bigsi = if high_mem_load {
                bigsi::read_bigsi_highmem(index)
            } else {
                bigsi::read_bigsi(index)
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
            if fq[0].ends_with(".gz") {
                if fq.len() > 1 {
                    colorid::read_id_mt_pe::per_read_stream_pe(
                        fq,
                        &bigsi.map,
                        &bigsi.colors,
                        &bigsi.n_ref_kmers,
                        bigsi.bloom_size,
                        bigsi.num_hash,
                        bigsi.k_size,
                        0,
                        threads,
                        down_sample,
                        fp_correct,
                        batch,
                        prefix,
                        quality,
                    )
                } else {
                    colorid::read_id_mt_pe::per_read_stream_se(
                        fq,
                        &bigsi.map,
                        &bigsi.colors,
                        &bigsi.n_ref_kmers,
                        bigsi.bloom_size,
                        bigsi.num_hash,
                        bigsi.k_size,
                        0,
                        threads,
                        down_sample,
                        fp_correct,
                        batch,
                        prefix,
                        quality,
                    )
                };
            } else {
                colorid::read_id_mt_pe::stream_fasta(
                    fq,
                    &bigsi.map,
                    &bigsi.colors,
                    &bigsi.n_ref_kmers,
                    bigsi.bloom_size,
                    bigsi.num_hash,
                    bigsi.k_size,
                    0,
                    threads,
                    down_sample,
                    fp_correct,
                    batch,
                    prefix,
                );
            }
        }
        colorid::reports::read_counts_five_fields(prefix.to_owned() + "_reads.txt", prefix);
    }
    if let Some(matches) = matches.subcommand_matches("read_filter") {
        let classification = matches.value_of("classification").unwrap();
        let files: Vec<_> = matches.values_of("files").unwrap().collect();
        let taxon = matches.value_of("taxon").unwrap();
        let prefix = matches.value_of("prefix").unwrap();
        let exclude = matches.is_present("exclude");
        let map = colorid::read_filter::tab_to_map(classification.to_string(), taxon);
        if files.len() == 1 {
            colorid::read_filter::read_filter_se(map, files, taxon, prefix, exclude);
        } else {
            colorid::read_filter::read_filter_pe(map, files, taxon, prefix, exclude);
        }
    }
}
