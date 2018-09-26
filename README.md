[![Build Status](https://travis-ci.org/hcdenbakker/bigs_id.svg?branch=master)](https://travis-ci.org/hcdenbakker/bigs_id)

# colorid

An experiment with writing code in Rust and the BIGSI data-structure (preprint: https://www.biorxiv.org/content/early/2017/12/18/234955) and Phelim Bradley's implementation (https://github.com/Phelimb/BIGSI).

## Download a precompiled binary
Here: https://github.com/hcdenbakker/colorid/releases

## or build it yourself with Rust!

Install Rust on your system (https://www.rust-lang.org/en-US/install.html)

Clone this repository:

```git clone https://github.com/hcdenbakker/colorid.git```

Get into the bigs_id directory:
```cd colorid```

And build your binary:

```cargo build --release```

The binary can now be found in the `/target/release` directory within the `colorid` directory. Add the binary to your path for easy access.

## Usage
```
FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

SUBCOMMANDS:
    build          builds a bigsi
    help           Prints this message or the help of the given subcommand(s)
    info           dumps index parameters and accessions
    read_filter    filters reads
    read_id        id's reads
    search         does a bigsi search on one or more fasta/fastq.gz file
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

``` ./target/release/colorid build -r ref_file_example.txt -b test -k 31 -s 50000000 -n 4```

Note! These parameters work well for single isolate, mixed samples with a few species. For complex metagenomic samples the BIGSI parameters need to be adjusted to adjust the false positive rate.

### 2. Search

``` ./target/release/colorid search -b test.bxi -q SRR4098796_1.fastq.gz -r SRR4098796_2.fastq.gz ```

### 3. results
With the default settings `bigs_id` will report reference sequences that share >35% of their k-mers with the query (more about this threshold to follow later). Here is the output of a search with SRA accession SRR4098796 (L. monocytogenes lineage I) as query:
```
SRR4098796_1.fastq.gz	3076072	Listeria_monocytogenes_F2365	0.87	134.25	126	475266
SRR4098796_1.fastq.gz	3076072	Listeria_monocytogenes_SRR2167842	0.40	128.25	122	7831
```
In the first column we find the query, the second column shows the number of k-mers in the query, the third column displays the reference sequence, the fourth column the proportion of kmers in the reference shared with the query, the fifth column displays the average coverage based on k-mers that were uniquely matched with this reference, the sixth the modus of the coverage based on uniquely matched k-mers and the last column the number of uniquely matched k-mers.

## Interrogate unassembled genome data for the presence of specific genes:

This is the main use case described in the biorxiv pre-print. Implemented are two algorithms, the `-g` algorithm, which allows for imperfect matches, and the fast `-s` algorithm, which reports only perfect matches.

In this case I assume we have couple of hundred unassembled Listeria monocytogenes genomes, and we want to check for the presence a specific gene in our sequenced genomes.

First we need to make a 'reference file' for the files we want to query. I like to use a simple bash loop to do this:
```
 for f in *_L001_R1_001.fastq.gz; do echo ${f%_L001_R1_001.fastq.gz}$'\t'$f >> example.txt; done
```

Next we create an index of the unassembled genome data. Given the fact our query size will be relatively small, we can use much stringent parameters for the index:

``` ./target/release/colorid build -r example.txt -b 30M_2H_K21.bxi -k 21 -s 30000000 -n 2```

This will perform a single threaded build of the index. If you want to speed up the build, you can use the `-t` flag to run it in multithreaded mode: 

``` ./target/release/colorid build -r example.txt -b 30M_2H_K21.bxi -k 21 -s 30000000 -n 2 -t 24```
This command will perform a build in 24 threads.

Now it is time for the search!

``` ./target/release/colorid search -b 30M_2H_K21.bxi -q geneX.fasta -g ```

When we use this subcommand we get results presented as follows:
```
geneX.fasta	ex1	1143	0.964
geneX.fasta	ex2	1143	1.000
geneX.fasta	ex3	1143	1.000
geneX.fasta	ex4	1143	0.982
geneX.fasta	ex5	1143	0.982
```

The first column gives us the query name, the second accesion in the index with a hit, the third the number of kmers in the query, and the fourth the proportion of k-mers of our query that were found in this accession. The `-s` flag will present the results in a similar way, only presenting perfect matches.  

## Classifying reads with read_id

This subcommand uses a simple majority-rule algorithm to classify reads, and the results can be used with the read_filter subcommand to either create a read file for a specific taxon, or filter a specific taxon from a read file. Write up to follow soon! 

## Acknowledgements
Lee Katz (https://github.com/lskatz) for his help with setting up Travis CI! 
