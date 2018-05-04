[![Build Status](https://travis-ci.org/hcdenbakker/bigs_id.svg?branch=master)](https://travis-ci.org/hcdenbakker/bigs_id)

# bigs_id

An experiment with writing code in Rust and the BIGSI data-structure (preprint: https://www.biorxiv.org/content/early/2017/12/18/234955) and Phelim Bradley's implementation (https://github.com/Phelimb/BIGSI).

## Download a precompiled binary
Here: https://github.com/hcdenbakker/bigs_id/releases

## or build it yourself with Rust!

Install Rust on your system (https://www.rust-lang.org/en-US/install.html)

Clone this repository:

```git clone https://github.com/hcdenbakker/bigs_id.git```

Get into the bigs_id directory:
```cd bigs_id```

And build your binary:

```cargo build --release```

The binary can now be found in the `/target/release` directory within the `bigs_id` directory. Add the binary to your path for easy access.

## Usage
```
FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

SUBCOMMANDS:
    build      builds a bigsi
    help       Prints this message or the help of the given subcommand(s)
    info       dumps index parameters and accessions
    read_id    id's reads
    search     does a bigsi search on one or more fasta/fastq.gz files
```

# Example:

## Create index

``` ./target/release/bigs_id build -r ref_file_example.txt -b test.bxi -k 31 --s 50000000 -n 4```

Note! These parameters work well for single isolate, mixed samples with a few species. For complex metagenomic samples the BIGSI parameters need to be adjusted to adjust the false positive rate.

## Search

``` ./target/release/bigs_id search -b test.bxi -q your_query.fastq.gz ```

Takes gzipped fastq or fasta files

More examples of uses to follow later.

Enjoy!

## Acknowledgements
Lee Katz (https://github.com/lskatz) for his help with setting up Travis CI! 
