conda activate iss

# to get some non uniformity
iss generate --genomes Sa2_SA222.fasta --cpus  8 --model novaseq --compress -n 75000 -o Sa2int -C halfnormal

# generate some background
iss generate --genomes C222.fasta  --cpus  8 --model novaseq --compress -n 75000 -o C222 -C exponential


 cat C222_R1.fastq.gz Sa2int_R1.fastq.gz > combo_R1.fastq.gz
 cat C222_R2.fastq.gz Sa2int_R2.fastq.gz > combo_R2.fastq.gz

