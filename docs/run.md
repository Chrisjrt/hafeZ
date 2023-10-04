# Running hafeZ

Once `hafeZ database` has been run, hafeZ requires an input genome FASTA file and sequencing read FASTQ files(s). An output directory can be specified using -o. Otherwise, an `output/` directory will be created in your current working directory.

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


