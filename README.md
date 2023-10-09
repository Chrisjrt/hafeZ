<p align="center">
  <img src="hafeZ_logo.jpg" width="40%">
</p>

# hafeZ
A tool for identifying inducible prophage elements through read mapping

"*I caught the happy virus last night when I was out singing beneath the stars.*"
-Hafez

# Quick Start

The easiest way to install hafeZ is via conda:

`conda install -c bioconda hafez`

Followed by database installation using:

`hafeZ database`

And finally running hafeZ `short` or `long`:

`hafeZ short -g genome.fasta -1 read_1.fastq.gz -2 read_2.fastq.gz -o output_folder -t <threads>`
`hafeZ long -g genome.fasta -1 read_1.fastq.gz -2 read_2.fastq.gz -o output_folder -t <threads>`


# Table of Contents

- [hafeZ](#hafez)
- [Quick Start](#quick-start)
- [Table of Contents](#table-of-contents)
- [Installation](#installation)
  - [Bioconda](#bioconda)
  - [Source](#source)
    - [Python dependencies (will be installed with the pip command)](#python-dependencies-will-be-installed-with-the-pip-command)
    - [Other dependencies (need to be installed separately)](#other-dependencies-need-to-be-installed-separately)
- [Quick start](#quick-start-1)
  - [Help](#help)
  - [Database download setup](#database-download-setup)
  - [Short Reads](#short-reads)
  - [Long Reads](#long-reads)
  - [Test dataset](#test-dataset)
  - [Outputs](#outputs)
- [Citation](#citation)

# Installation

## Bioconda

```
mamba create -n hafeZ -c conda-forge -c bioconda -c defaults hafez
```

## Source

```
git clone "https://github.com/Chrisjrt/hafeZ.git"
cd hafeZ
pip install -e .
hafeZ -h
```

If installing from source, hafeZ requires the following dependencies to also be installed:

### Python dependencies (will be installed with the pip command)

- pandas
- numpy
- matplotlib
- scipy
- Biopython
- [pyrodigal](https://github.com/althonos/pyrodigal)
- pysam
- [pyhmmer](https://github.com/althonos/pyhmmer)
- loguru
- alive-progress
- requests

### Other dependencies (need to be installed separately)

- [minimap2](https://github.com/lh3/minimap2)
- [samtools](https://github.com/samtools)
- [mosdepth](https://github.com/brentp/mosdepth)
- [rasusa](https://github.com/mbhall88/rasusa)
- [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi)

# Quick start

## Help

To access the help menu use the `-h` option:

```
Usage: hafeZ [OPTIONS] COMMAND [ARGS]...

Options:
  -h, --help     Show this message and exit.
  -V, --version  Show the version and exit.

Commands:
  citation  Print the citation for hafeZ
  database  Downloads and installs hafeZ database
  long      Runs hafeZ with ONT long reads
  short     Runs hafeZ with paired end short reads
```


## Database download setup

hafeZ uses the PHROGs HMM database. It can be downloaded and installed to the default location using

```
hafeZ database
```

If you want to specify a specific directory to install the database, use `-d`

```
hafeZ database -d <database directory>
```

## Short Reads

`hafeZ short` takes paired end short reads in both `.fastq` and `.fastq.gz` format. 

```
hafeZ short -g genome.fasta -1 read_1.fastq.gz -2 read_2.fastq.gz -o output_folder 
```

A full list of parameters can be found below.

```
Usage: hafeZ short [OPTIONS]

  Runs hafeZ with paired end short reads

Options:
  -h, --help                    Show this message and exit.
  -V, --version                 Show the version and exit.
  -1, --reads1 PATH             Path to R1 reads file in FASTQ or FASTQ gzip
                                format.  [required]
  -2, --reads2 PATH             Path to R2 reads file in FASTQ or FASTQ gzip
                                format  [required]
  -g, --genome PATH             Path to genome assembly in FASTA format
                                [required]
  -o, --output PATH             Output directory path  [default: output_hafeZ]
  -t, --threads INTEGER         Number of threads.  [default: 1]
  -d, --database TEXT           Path to the hafeZ database directory.
                                [default: database]
  -f, --force                   Force overwrites the output directory
  -b, --bin_size INTEGER        Bin size in bp to use for coverage depth
                                smoothing. Must be an odd number.  [default:
                                3001]
  -c, --cutoff FLOAT            Z-score cutoff for initially detecting RoIs.
                                [default: 3.5]
  -w, --width INTEGER           Minimum width (bp) of RoI that passes Z-score
                                cutoff.  [default: 4000]
  -m, --min_orfs INTEGER        Minimum number of ORFs needed in an RoI.
                                [default: 5]
  -p, --phrog_fract FLOAT       Minimum fraction number of ORFs needed in an
                                RoI with PHROG hit.  [default: 0.1]
  -z, --median_z_cutoff FLOAT   Median Z-score for an roi to be retained.
                                [default: 3.5]
  -S, --sub_sample              Randomly sub-sample reads to adjust overall
                                mapping coverage of genome. N.B. Use
                                -C/--coverage to adjust coverage desired.
  -C, --coverage FLOAT          Desired coverage of genome to be used for
                                subsampling reads (default = 100.0). N.B. Must
                                be used with -S\--sub_sample flag.  [default:
                                100.0]
  -N, --no_extra                Turn off extra accuracy checks using clipped
                                sequences. Might give same results, might give
                                extra rois.
  -M, --memory_limit TEXT       Upper bound per thread memory limit for
                                samtools, suffix K/M/G recognized (default =
                                768M).  [default: 768M]
  -Z, --all_zscores             Make graphs of all contig Z-scores even if no
                                roi found (useful for manual inspection). N.B.
                                This will make graphs for each contig, so if
                                you have 100 contigs you will get 100 graphs.
  --min_contig_len INTEGER      Minimum contig length that hafeZ will consider
                                as potentially harbouring a prophage.
                                [default: 10000]
  --join_window INTEGER         Minimum window within which 2 ROIs will be
                                merged.  [default: 10000]
  --evalue FLOAT                Evalue threshold for significant PyHMMER hits.
                                [default: 0.001]
  -e, --expect_mad_zero         allow MAD == 0 to exit without non-zero exit
                                code. Will also cause coverage plots for each
                                contig to be output to help with debugging.
                                Useful for uninduced lysates.
  -k, --keep_threshold INTEGER  Threshold for number of best soft clip
                                combinations to keep for each roi.  [default:
                                50]
```



## Long Reads

`hafeZ long` takes single end (ONT) long reads in both `.fastq` and `.fastq.gz` format. 

```
hafeZ long -g genome.fasta -l longreads.fastq.gz -o output_folder 
```

A full list of parameters can be found below.

```
Usage: hafeZ long [OPTIONS]

  Runs hafeZ with ONT long reads

Options:
  -h, --help                   Show this message and exit.
  -V, --version                Show the version and exit.
  -l, --longreads PATH         Path to longreads file in FASTQ or FASTQ gzip
                               format.  [required]
  -g, --genome PATH            Path to genome assembly in FASTA format
                               [required]
  -o, --output PATH            Output directory path  [default: output_hafeZ]
  -t, --threads INTEGER        Number of threads.  [default: 1]
  -d, --database TEXT          Path to the hafeZ database directory.
                               [default: database]
  -f, --force                  Force overwrites the output directory
  -b, --bin_size INTEGER       Bin size in bp to use for coverage depth
                               smoothing. Must be an odd number.  [default:
                               3001]
  -c, --cutoff FLOAT           Z-score cutoff for initially detecting RoIs.
                               [default: 3.5]
  -w, --width INTEGER          Minimum width (bp) of RoI that passes Z-score
                               cutoff.  [default: 4000]
  -m, --min_orfs INTEGER       Minimum number of ORFs needed in an RoI.
                               [default: 5]
  -p, --phrog_fract FLOAT      Minimum fraction number of ORFs needed in an
                               RoI with PHROG hit.  [default: 0.1]
  -z, --median_z_cutoff FLOAT  Median Z-score for an roi to be retained.
                               [default: 3.5]
  -S, --sub_sample             Randomly sub-sample reads to adjust overall
                               mapping coverage of genome. N.B. Use
                               -C/--coverage to adjust coverage desired.
  -C, --coverage FLOAT         Desired coverage of genome to be used for
                               subsampling reads (default = 100.0). N.B. Must
                               be used with -S\--sub_sample flag.  [default:
                               100.0]
  -N, --no_extra               Turn off extra accuracy checks using clipped
                               sequences. Might give same results, might give
                               extra rois.
  -M, --memory_limit TEXT      Upper bound per thread memory limit for
                               samtools, suffix K/M/G recognized (default =
                               768M).  [default: 768M]
  -Z, --all_zscores            Make graphs of all contig Z-scores even if no
                               roi found (useful for manual inspection). N.B.
                               This will make graphs for each contig, so if
                               you have 100 contigs you will get 100 graphs.
  --min_contig_len INTEGER     Minimum contig length that hafeZ will consider
                               as potentially harbouring a prophage.
                               [default: 10000]
  --join_window INTEGER        Minimum window within which 2 ROIs will be
                               merged.  [default: 10000]
  --evalue FLOAT               Evalue threshold for significant PyHMMER hits.
                               [default: 0.001]
  -e, --expect_mad_zero        allow MAD == 0 to exit without non-zero exit
                               code. Will also cause coverage plots for each
                               contig to be output to help with debugging.
                               Useful for uninduced lysates.
  --min_reads INTEGER          Minimum number of longreads that need to map to
                               coordinates on both ends of an ROI to be
                               classified as induced.  [default: 5]
```



## Test dataset

You can run hafeZ on a test datasets found in the `tests` directory based on simulated data of _S. aureus_ C222 isolate (from [here](https://www.biorxiv.org/content/10.1101/2023.03.28.534496v1)).

Alternatively, you can run it on a full dataset:

```
mkdir test
wget -P test/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR455/005/ERR4552545/ERR4552545_1.fastq.gz
wget -P test/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR455/005/ERR4552545/ERR4552545_2.fastq.gz
wget -P test/ https://www.ebi.ac.uk/ena/browser/api/fasta/CP015406.2?download=true
mv test/CP015406.2?download=true test/CP015406.2.fasta
# assuming you have the database installed and the environment set up
hafeZ short -r1 test/ERR4552545_1.fastq.gz -r2 test/ERR4552545_2.fastq.gz -o test/output -g test/CP015406.2.fasta -t 8 
```

## Outputs

If a putative active prophage is found hafeZ produces six main ouputs:

- `hafeZ_all_roi_seqs.fasta` = file containing the DNA sequences of all the regions of interest identified
- `hafeZ_summary_all_rois.tsv`  = file containing a summary of info related to each region of interest
- `hafeZ_pyhmmer_hits.tsv` = file containing a list of all region if interest orfs and a description of the pvogs/phrogs they hit
- `hafeZ_prophages_for_xxx.png` = image of zscores per base within contigs where a region of interest was identified with the region highlights (one file per contig containing a region of interest)
- `hafeZ_orfs_aa_XXX.fasta` = fasta file containing amino acid sequence of each orf within the roi
- `hafeZ_orfs_dna_XXX.fasta` = fasta file containing the dna sequence of each orf within the roi

N.B. if the `-Z` flag is used an additional input, `zscores_for_contigXXX.png`, will also be generated which shows the Z-scores of each contig examined (i.e. if input genome contains 100 contigs there will be 100 zscore .png files output. This can be useful if the user wants to manuallly inspect for any potential rois that may be missed under default paramaters. )

If no putative active prophages are found hafeZ will output only an empty `hafeZ_summary_all_rois.tsv` file. 

Additionally, the `hafez__.log` file will contain the logging output, and `logs` directory will contain the logs for all the subprocesses in hafeZ (minimap2, samtools, blast etc).


# Citation

If you publish results from hafeZ please cite the following:

hafeZ: Active prophage identification through read mapping (bioRxiv)
https://doi.org/10.1101/2021.07.21.453177

https://github.com/Chrisjrt/hafeZ
