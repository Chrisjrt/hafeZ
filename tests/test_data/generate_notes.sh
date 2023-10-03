conda activate iss

# to get some non uniformity
iss generate --genomes Sa2_SA222.fasta --cpus  8 --model novaseq --compress -n 75000 -o Sa2int -C halfnormal

# generate some background
iss generate --genomes C222.fasta  --cpus  8 --model novaseq --compress -n 75000 -o C222 -C exponential


 cat C222_R1.fastq.gz Sa2int_R1.fastq.gz > combo_R1.fastq.gz
 cat C222_R2.fastq.gz Sa2int_R2.fastq.gz > combo_R2.fastq.gz




# generate some long reads

conda activate badread

badread simulate --reference Sa2_SA222.fasta  --quantity 50x  --error_model nanopore2023 --qscore_model nanopore2023 --seed 43 \
         --length 10000,10000 | gzip > Sa2int.fastq.gz


badread simulate --reference C222.fasta  --quantity 0.5x  --error_model nanopore2023 --qscore_model nanopore2023 --seed 43 \
         --length 10000,10000 | gzip > C222.fastq.gz

# tests

hafeZ short -g tests/test_data/short/C222.fasta -1 tests/test_data/short/combo_R1.fastq.gz -2 tests/test_data/short/combo_R2.fastq.gz -o test_out -t 8 -f -C 15 -S -d ../pharokka_v1.4.0_databases --all_zscores

hafeZ long -g tests/test_data/long/C222.fasta -l tests/test_data/long/combo.fastq.gz -o test_out_long -t 8 -f -C 3.6 -S -d ../pharokka_v1.4.0_databases --all_zscores