<p align="center">
  <img src="hafeZ_logo.jpg" width="40%">
</p>

# hafeZ
A tool for identifying inducible prophage elements through read mapping

"*I caught the happy virus last night when I was out singing beneath the stars.*"
-Hafez

# Installation

## Bioconda

```
conda env create --name hafeZ --file=environment.yml
```

## Source

If installing from source, hafeZ requires the following dependencies to also be installed:

### Python dependencies

- pandas
- numpy
- matplotlib
- scipy
- Biopython
- pyrodigal
- pysam

### Other dependencies

- minimap2
- samtools
- mosdepth
- hmmer3
- blast

# Quick start

## Help

To access the help menu use the `-h` option:

```
hafeZ.py -h
```


## Initial setup

As hafeZ uses the pVOGs database this must first be retrieved and formatted before use. This can be done using the following command:

```
hafeZ.py -G hafeZ_db/
```

## illumina reads

hafeZ accepts illumina reads in both .fastq and .fastq.gz format. To use hafeZ with illumina reads:

```
hafeZ.py -f assembly.fasta -r1 read_1.fastq.gz -r2 read_2.fastq.gz -o output_folder -D hageZ_db
```

## Test dataset

A test dataset can be obtained and ran using the following:

```
mkdir test
wget -P test/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR455/005/ERR4552545/ERR4552545_1.fastq.gz
wget -P test/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR455/005/ERR4552545/ERR4552545_2.fastq.gz
wget -P test/ https://www.ebi.ac.uk/ena/browser/api/fasta/CP015406.2?download=true
mv test/CP015406.2?download=true test/CP015406.2.fasta
./hafeZ.py -r1 test/ERR4552545_1.fastq.gz -r2 test/ERR4552545_2.fastq.gz -o test/output -O -f test/CP015406.2.fasta -t 8 -D hafeZ_db/ -Z
```

## Outputs

hafeZ produces five main ouputs:

- hafeZ_all_roi_seqs.fasta = file containing the DNA sequences of all the regions of interest identified
- hafeZ_all_rois.tsv = file containing a summary of info related to each region of interest
- hafeZ_prophage_for_xxx.png = image of zscores per base within contigs where a region of interest was identified with the region highlights (one file per contig containing a region of interest)
- hafeZ_orfs_aa_XXX.faa = fasta file containing amino acid sequence of each orf within the roi
- hafeZ_orfs_dna_XXX.fasta = fasta file containing the dna sequence of each orf within the roi

N.B. if the -Z option is used an additional input, zscores_for_contigXXX.png, will also be generated which shows the Z-scores of each contig examined (i.e. if input genome contains 100 contigs there will be 100 zscore .png files output. This can be useful if the user wants to manuallly inspect for any potential rois that may be missed under default paramaters. )

## Caveat

hafeZ is currently only optimised to use paired end illumina reads as inputs. Future updates will allow use of single end illumina reads, nanopore reads, and pacbio reads, but these have not optimised yet.

# Citation

If you publish results from hafeZ please cite the following:

doi.org/some.preprint

https://github.com/Chrisjrt/hafeZ
