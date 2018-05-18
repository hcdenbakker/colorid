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

# Example uses:

## Infer k-mer similarity between query and reference sequences:

Infer k-mer similarity between a query fasta/fastq.gz file and an index of reference sequences with the `search` subcommand:
```
-b, --bigsi=[FILE] 'Sets the name of the index file for search'
-q, --query      'one or more fastq.gz formatted files to be queried'
-f, --filter     'Sets threshold to filter k-mers by frequency'
-p, --p_shared        'minimum proportion of kmers shared with reference'
-c, --compressed 'If set to 'true', will assume compressed index (default: false)'
-g, If set('-g'), the proportion of kmers from the query matching the entries in the index will be reported'
-s, If set('-s'), the 'speedy' 'perfect match' alforithm will be performed'
```

### 1. Create index

``` ./target/release/bigs_id build -r ref_file_example.txt -b test.bxi -k 31 -s 50000000 -n 4```

Note! These parameters work well for single isolate, mixed samples with a few species. For complex metagenomic samples the BIGSI parameters need to be adjusted to adjust the false positive rate.

### 2. Search

``` ./target/release/bigs_id search -b test.bxi -q your_query.fastq.gz ```

### 3. results
With the default settings `bigs_id` will report reference sequences that share 35% of their k-mers with the query (more about this threshold to follow later). Here is the output of a search with SRA accession SRR4098796 (L. monocytogenes lineage I) as query:
```
Listeria_monocytogenes_F2365: 0.87 68.43 64 664676
Listeria_monocytogenes_SRR2167842: 0.40 62.65 2 11416
```
In the first column we find the reference sequence, the second column shows the proportion of k-mers in the reference shared with the query, the third column displays the average coverage based on k-mers that were uniquely matched with this reference, the fourth the modus of the coverage based on uniquely matched k-mers and the last column the number of uniquely matched k-mers.

## Interrogate unassembled genome data for the presence of specific genes:

This is the main use case described in the biorxiv pre-print. Implemented are two algorithms, the `-g` algorithm, which allows for imperfect matches, and the fast `-s` algorithm, which reports only perfect matches.

In this case I assume we have couple of hundred unassembled Listeria monocytogenes genomes, and we want to check for the presence a cadmium resistance gene (e.g., CadB) in our sequenced genomes.

First we need to make a 'reference file' for the files we want to query. I like to use a simple bash loop to do this:
```
 for f in *_L001_R1_001.fastq.gz; do echo ${f%_L001_R1_001.fastq.gz}$'\t'$f >> example.txt; done
```

Next we create an index of the unassembled genome data. Given the fact our query size will be relatively small, we can use much stringent parameters for the index:

``` ./target/release/bigs_id build -r example.txt -b 30M_2H_K21.bxi -k 21 -s 30000000 -n 2```

This will perform a single threaded build of the index. If you want to speed up the build, you can use the `-t` flag to run it in multithreaded mode: 

``` ./target/release/bigs_id build -r example.txt -b 30M_2H_K21.bxi -k 21 -s 30000000 -n 2 -t 24```
This command will perform a build in 24 threads.

Now it is time for the search!

``` ./target/release/bigs_id search -b 30M_2H_K21.bxi -q pLM33_CadB.fasta -g ```

When we use this subcommand we get results presented as follows:
```ex1: 1.00
ex2: 0.99
ex3: 0.97
ex4: 0.74
ex5: 1.00
```

The first column gives us the accesion in the index, the second the proportion of k-mers of our query that were found in this accession. The `-s' results just consist of one column, the accessions for which a perfect match for the query was found.

## Classifying reads with read_id

This subcommand uses a simple majority-rule algorithm to classify reads. Write up to follow soon! 

## Acknowledgements
Lee Katz (https://github.com/lskatz) for his help with setting up Travis CI! 
