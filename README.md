# bigs_id

## build

```cargo build --release```

## Usage
```USAGE:
    bigs_id [SUBCOMMAND]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

SUBCOMMANDS:
    batch_search    does a bigsi search on a bunch of fastq.gz files
    build           builds a bigsi
    help            Prints this message or the help of the given subcommand(s)
    search          does a bigsi search
```

## Create index

``` ./target/release/bigs_id build -r ref_file_example.txt -b test.bxi -k 31 --bloom 50000000 -n 4```

Note! These parameters work well for single isolate, mixed samples with a few species. For complex metagenomic samples the BIGSI parameters need to be adjusted to adjust the false positive rate.

## Search

``` ./target/release/bigs_id search -b test.bxi -q your_query.fastq.gz ```

Takes gzipped fastq or fasta files

Enjoy
