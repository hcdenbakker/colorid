./target/release/bigs_id build -s 750000 -n 4 -k 27 -b ./test_data/phage.bxi -r ./test_data/ref_file.txt

./target/release/bigs_id search -b ./test_data/phage.bxi -q ./test_data/SRR548019.fastq.gz -f 1
