#!/usr/bin/env python3
"""hafeZ"""

import os
import sys
from pathlib import Path

import click
import numpy as np
import pandas as pd
from loguru import logger

from hafeZ.utils.db import check_db_installation
from hafeZ.utils.exit import exit_error_gracefully, exit_success
from hafeZ.utils.extra_process_roi import quick_filter
from hafeZ.utils.mapping import get_bam, get_cov, minimap_long, minimap_paired
from hafeZ.utils.mapping_calcs import get_ZScores, smooth_depths
from hafeZ.utils.orfs import (
    calc_phrogs_frac,
    extract_roi_orfs,
    get_orfs,
    get_roi_sequences,
    run_pyhmmer,
)
from hafeZ.utils.post_processing import (
    get_att,
    get_names,
    output_all_phrogs,
    output_contig_Z,
    output_prophage_graphs,
    output_roi_orfs,
    output_roi_seqs,
    output_roi_table,
)
from hafeZ.utils.premap import Premap
from hafeZ.utils.process_rois import Haf
from hafeZ.utils.util import begin_hafeZ, get_version, print_citation
from hafeZ.utils.validation import (
    check_dependencies,
    check_memory_limit,
    instantiate_dirs,
    validate_fasta,
)

"""
some code adapted from tbpore https://github.com/mbhall88/tbpore
"""

log_fmt = (
    "[<green>{time:YYYY-MM-DD HH:mm:ss}</green>] <level>{level: <8}</level> | "
    "<level>{message}</level>"
)


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option(
            "-g",
            "--genome",
            help="Path to genome assembly in FASTA format",
            type=click.Path(),
            required=True,
        ),
        click.option(
            "-o",
            "--output",
            default="output_hafeZ",
            show_default=True,
            type=click.Path(),
            help="Output directory path",
        ),
        click.option(
            "-t",
            "--threads",
            help="Number of threads.",
            default=1,
            show_default=True,
        ),
        click.option(
            "-d",
            "--database",
            help="Path to the hafeZ database directory.",
            default="database",
            show_default=True,
        ),
        click.option(
            "-f", "--force", is_flag=True, help="Force overwrites the output directory"
        ),
        click.option(
            "-b",
            "--bin_size",
            type=int,
            default=3001,
            help="Bin size in bp to use for coverage depth smoothing. Must be an odd number.",
            show_default=True,
        ),
        click.option(
            "-c",
            "--cutoff",
            type=float,
            default=3.5,
            help="Z-score cutoff for initially detecting RoIs.",
            show_default=True,
        ),
        click.option(
            "-w",
            "--width",
            type=int,
            default=4000,
            help="Minimum width (bp) of RoI that passes Z-score cutoff.",
            show_default=True,
        ),
        click.option(
            "-m",
            "--min_orfs",
            type=int,
            default=5,
            help="Minimum number of ORFs needed in an RoI.",
            show_default=True,
        ),
        click.option(
            "-p",
            "--phrog_fract",
            type=float,
            default=0.1,
            help="Minimum fraction number of ORFs needed in an RoI with PHROG hit.",
            show_default=True,
        ),
        click.option(
            "-z",
            "--median_z_cutoff",
            type=float,
            default=3.5,
            help="Median Z-score for an roi to be retained.",
            show_default=True,
        ),
        click.option(
            "-S",
            "--sub_sample",
            is_flag=True,
            help="Randomly sub-sample reads to adjust overall mapping coverage of genome. N.B. Use -C/--coverage to adjust coverage desired.",
        ),
        click.option(
            "-C",
            "--coverage",
            type=float,
            default=100.0,
            help="Desired coverage of genome to be used for subsampling reads (default = 100.0). N.B. Must be used with -S\--sub_sample flag.",
            show_default=True,
        ),
        click.option(
            "-N",
            "--no_extra",
            is_flag=True,
            help="Turn off extra accuracy checks using clipped sequences. Might give same results, might give extra rois.",
        ),
        click.option(
            "-M",
            "--memory_limit",
            type=str,
            default="768M",
            help="Upper bound per thread memory limit for samtools, suffix K/M/G recognized (default = 768M).",
            show_default=True,
        ),
        click.option(
            "-Z",
            "--all_zscores",
            is_flag=True,
            help="Make graphs of all contig Z-scores even if no roi found (useful for manual inspection). N.B. This will make graphs for each contig, so if you have 100 contigs you will get 100 graphs.",
        ),
        click.option(
            "--min_contig_len",
            type=int,
            default=10000,
            help="Minimum contig length that hafeZ will consider as potentially harbouring a prophage.",
            show_default=True,
        ),
        click.option(
            "--join_window",
            type=int,
            default=10000,
            help="Minimum window within which 2 ROIs will be merged.",
            show_default=True,
        ),
        click.option(
            "--evalue",
            type=float,
            default=0.001,
            help="Evalue threshold for significant PyHMMER hits.",
            show_default=True,
        ),
        click.option(
            "-e",
            "--expect_mad_zero",
            is_flag=True,
            help="allow MAD == 0 to exit without non-zero exit code. Will also cause coverage plots for each contig to be output to help with debugging. Useful for uninduced lysates.",
        ),
    ]
    for option in reversed(options):
        func = option(func)
    return func


@click.group()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
def main_cli():
    1 + 1


"""
short command
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@click.option(
    "-1",
    "--reads1",
    help="Path to R1 reads file in FASTQ or FASTQ gzip format.",
    type=click.Path(),
    required=True,
)
@click.option(
    "-2",
    "--reads2",
    help="Path to R2 reads file in FASTQ or FASTQ gzip format",
    type=click.Path(),
    required=True,
)
@common_options
@click.option(
    "-k",
    "--keep_threshold",
    type=int,
    default=50,
    help="Threshold for number of best soft clip combinations to keep for each roi.",
    show_default=True,
)
def short(
    ctx,
    genome,
    reads1,
    reads2,
    output,
    threads,
    database,
    force,
    bin_size,
    cutoff,
    width,
    min_orfs,
    phrog_fract,
    median_z_cutoff,
    keep_threshold,
    sub_sample,
    coverage,
    no_extra,
    memory_limit,
    all_zscores,
    min_contig_len,
    join_window,
    expect_mad_zero,
    evalue,
    **kwargs,
):
    """Runs hafeZ with paired end short reads"""

    # validates the directory  (need to before I start hafez or else no log file is written)
    instantiate_dirs(output, force)

    logdir: Path = Path(output) / "logs"

    output: Path = Path(output)

    params = {
        "--genome": genome,
        "--reads1": reads1,
        "--reads2": reads2,
        "--output": output,
        "--threads": threads,
        "--database": database,
        "--force": force,
        "--bin_size": bin_size,
        "--cutoff": cutoff,
        "--width": width,
        "--min_orfs": min_orfs,
        "--phrog_fract": phrog_fract,
        "--median_z_cutoff": median_z_cutoff,
        "--keep_threshold": keep_threshold,
        "--sub_sample": sub_sample,
        "--coverage": coverage,
        "--no_extra": no_extra,
        "--memory_limit": memory_limit,
        "--all_zscores": all_zscores,
        "--min_contig_len": min_contig_len,
        "--expect_mad_zero": expect_mad_zero,
        "--evalue": evalue,
    }

    # initial logging and list all the params
    start_time = begin_hafeZ(params)

    # validates fasta
    validate_fasta(genome)

    # check dependencies
    check_dependencies()

    # checks database installation
    logger.info(f"Checking database installation at {database}.")
    check_db_installation(database, install_flag=False)

    logger.info(f"Checking memory limit of {memory_limit}.")
    check_memory_limit(memory_limit)

    # instanatiate the premap class with coverage too
    ## premap.py
    premap = Premap()
    premap.coverage = coverage
    premap.output = output
    premap.start_time = start_time

    # filter the input genome to get long enough contigs
    premap.process_fasta(genome, min_contig_len)

    # get genome length of contigs
    premap.get_genome_length()

    # estimate coverage
    logger.info("Estimating overall genome coverage.")
    premap.calculate_estimated_coverage_short(reads1, reads2)
    logger.info(f"Estimated genome coverage is {premap.estimated_coverage:.2f}x.")

    # subsample if coverage is too high and set the Paths to the reads
    reads1, reads2 = premap.subsample_short(reads1, reads2, logdir)

    #######################
    # mapping.py
    # no class needed
    ########################
    # run the mapping commands - no class as files written to temp
    logger.info("Mapping reads to filtered genome and sorting.")

    # map to the filtered FASTA
    seq_temp_out: Path = output / "temp_genome.fasta"
    minimap_paired(reads1, reads2, seq_temp_out, output, threads, logdir)
    get_bam(output, threads, memory_limit, logdir)
    logger.info("Calculating per base coverage.")
    get_cov(output, threads, logdir)

    #######################
    # mapping_calcs.py
    #######################
    # smooths depths
    ## note it fails if the induction is too clean!
    # i.e. if you just simulate the phage reads, so need a little bit of background

    logger.info("Smoothing signals.")
    depths, cov = smooth_depths(output, bin_size, threads)

    # get Z scores
    zscores, median, mad = get_ZScores(depths, output, start_time, cov, expect_mad_zero)

    ########################
    # ROI
    # process_roi.py
    # uses Haf class to store everything
    ########################

    # instanatiate the Haf class
    haf = Haf()
    haf.all_zscores = all_zscores
    haf.median = median
    haf.mad = mad
    haf.zscores = zscores
    haf.output = output
    haf.depths = depths
    haf.start_time = start_time

    logger.info("Finding Regions of Interest (ROIs).")
    haf.get_ROIs(cutoff, width)
    logger.info("Joining any ROIs that are close to each other.")
    # RoIs within this window will be joined
    haf.join_ROIs(join_window)

    # filter rois based on width
    logger.info("Removing any small ROIs.")
    haf.filter_ROIs(width)

    logger.info("Extracting depth information for each ROI.")
    # number of bases before and after the ROI
    neighbourhood = 15001
    logger.info("Parsing and processing the BAM file.")
    haf.bamfile_to_pandas(output, neighbourhood)

    # find ROIs that span an entire contig
    logger.info("Checking for ROIs near contig ends.")
    haf.find_contig_end_rois(neighbourhood)

    # find evidence of pairs of reads that are distant from each other (i.e. look to span deletion)
    logger.info("Finding distant read pairs.")
    haf.find_far_reads(
        width, neighbourhood
    )  ##### might want to make this a non exclusion stage and instead an informative stage (i.e. dont delete those that dont have, instead just make column value that says it )

    #### combine roi_df and end_df for next steps: ####
    if haf.end_roi_df is not None and haf.end_roi_df is not None:
        if len(haf.roi_df) > 0 and len(haf.end_roi_df) > 0:
            haf.roi_df = pd.concat([haf.roi_df, haf.end_roi_df])
            haf.roi_df = haf.roi_df.reset_index(drop=True)
        elif len(haf.roi_df) == 0 and len(haf.end_roi_df) > 0:
            haf.roi_df = haf.end_roi_df

    #### find soft clipped reads  ####
    if haf.end_roi_df is not None:
        if len(haf.roi_df) > 0:
            logger.info("Finding clipped reads.")
            haf.find_soft_clippings(width, threads, neighbourhood)

    #### refine rois at ends #
    if haf.end_roi_df is not None:
        if len(haf.roi_df) > 0:
            haf.find_contig_end_rois_again(neighbourhood)

    # calculate median Z-scores and filter rois

    if haf.end_roi_df is not None:
        if len(haf.roi_df) > 0:
            logger.info("Calculating median Z-score for each ROI.")
            #### calculate median Z-score for each roi ####
            haf.calc_roi_Z_medians(cov, median, mad, median_z_cutoff, threads)

    #### filter rois using median Z-score for each roi ####
    if haf.end_roi_df is not None:
        if len(haf.roi_df) > 0:
            logger.info("Filtering ROIs using the median Z-score.")
            haf.filter_Z_medians()

    ### filter to keep only best X rois ####
    if haf.end_roi_df is not None:
        if len(haf.roi_df) > 0:
            haf.keep_only_x_best(keep_threshold)

    #### check to see if any rois are circular ####
    if haf.end_roi_df is not None:
        if len(haf.roi_df) > 0:
            logger.info("Checking for circular ROIs (possible phage plasmids).")
            haf.check_for_circular_roi()

    # no extra
    if no_extra is False:
        ### collect the ends of the soft clipped reads found, map them, collect the results and process them ####

        if haf.end_roi_df is not None:
            if len(haf.roi_df) > 0:
                logger.info("Collecting clipped end reads.")
                haf.collecting_clipped_reads(output, threads)
                logger.info("Mapping clipped reads.")
                haf.map_clipped_reads(output, threads, genome, memory_limit, logdir)

                logger.info("Finding any ROIs that span an entire contig.")
                haf.seperate_rois_near_ends()

        if haf.end_roi_df is not None:
            if len(haf.roi_df) > 0:
                logger.info("Processing clipped read SAM file.")
                haf.process_clip_sam(output)

                if haf.clip_end_df.empty is True:  # only
                    exit_error_gracefully(
                        output, all_zscores, depths, median, mad, start_time
                    )
                    logger.error("No RoIs found at all with clipped reads.")
                else:
                    haf.get_clip_pos(output)

        # need to test this - my test data doesn't have it
        if haf.end_roi_df is not None:
            if len(haf.end_roi_df) > 0:
                haf.process_clip_sam_end_rois(output)
                haf.get_clip_pos_end_rois(premap.raw_seq_dict)

        if haf.end_roi_df is not None:
            if len(haf.end_roi_df) > 0:
                haf.filter_by_end_status()

        if haf.end_roi_df is not None:
            if len(haf.end_roi_df) > 0:
                haf.reformat_end_roi_tables(premap.raw_seq_dict)

    #### rejoin all databases together ####
    df_list = []
    for i in [haf.roi_df, haf.end_roi_df, haf.circular_df]:
        if i is not None:
            if len(i) > 0:
                df_list.append(i)
    if len(df_list) > 1:
        roi_df = pd.concat(df_list)
        roi_df = roi_df.reset_index(drop=True)
    elif len(df_list) == 1:
        roi_df = df_list[0]
    elif len(df_list) == 0:
        exit_error_gracefully(output, all_zscores, depths, median, mad, start_time)
        logger.error("No ROIs found.")

    # pick best rois if read checking has been done
    if no_extra is True:
        roi_df = quick_filter(roi_df)

    #### extract roi sequences ####
    for index, row in roi_df.iterrows():  #### move this elswhere
        if pd.isna(row["circular"]):
            roi_df.loc[index, "circular"] = False
        else:
            roi_df.loc[index, "med_z"] = np.nan

    ###########################
    # orfs.py
    ###########################

    ##### HMMs
    # roi sequences
    logger.info("Predicting CDS on all ROIs using Pyrodigal.")
    roi_dna = get_roi_sequences(roi_df, premap.filtered_seq_dict, output)

    #### call orfs in genome ####
    # orf df, AA and nt
    orf_df, orfs_aa, orfs_dna = get_orfs(premap.filtered_seq_dict, premap.multicontig)

    ### extract orfs for roi ###

    roi_df, roi_orf_aa, roi_orf_dna = extract_roi_orfs(
        orf_df,
        orfs_aa,
        roi_df,
        orfs_dna,
        output,
        min_orfs,
        all_zscores,
        depths,
        median,
        mad,
        start_time,
    )

    #### Screen all roi orfs vs phrogs db with pyhmmer ####

    evalue = 0.001
    evalue = 1.258e-60

    logger.info("Running PyHMMER on PHROGs.")
    best_phrog_results_dict = run_pyhmmer(database, output, threads, evalue=evalue)

    # calculate the percentage of orfs with PHROG hit for each roi
    logger.info("Calculating the proportion of ROI CDS with PHROGs hits.")
    roi_df = calc_phrogs_frac(
        best_phrog_results_dict,
        roi_df,
        phrog_fract,
        evalue,
        output,
        all_zscores,
        depths,
        median,
        mad,
        start_time,
    )

    #######################
    #### post_processing.py
    #######################

    # give rois new, final, names ahead of processing outputs

    roi_df = get_names(roi_df)
    logger.info("Finding possible attatchment sites.")
    roi_df = get_att(roi_df, premap.filtered_seq_dict, output, logdir)

    # save multifasta of roi genome seqs to file
    logger.info("Saving ROI sequences.")
    output_roi_seqs(roi_df, roi_dna, output)

    # output graph showing positions
    logger.info("Output plots.")
    output_prophage_graphs(roi_df, depths, output, median, mad)

    # output hmm table

    logger.info("Writing PyHMMER Output to file.")
    output_all_phrogs(roi_orf_aa, best_phrog_results_dict, output, database, evalue)

    #### output summary table of rois found ####
    logger.info("Writing ROI summary table.")
    output_roi_table(roi_df, output)

    # output roi orf aa and dna sequences
    logger.info("Writing ROI AA and DNA sequences.")
    output_roi_orfs(roi_orf_dna, roi_df, output, roi_orf_aa)

    # make Z-score graphs for each contig if -Z flag given
    if all_zscores is True:
        logger.info("Making Z-score graphs.")
        output_contig_Z(depths, output, median, mad)

    # exit hafeZ
    exit_success(output, start_time)


"""
long command
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@click.option(
    "-l",
    "--longreads",
    help="Path to longreads file in FASTQ or FASTQ gzip format.",
    type=click.Path(),
    required=True,
)
@common_options
@click.option(
    "--min_reads",
    help="Minimum number of longreads that need to map to coordinates on both ends of an ROI to be classified as induced.",
    type=int,
    default=5,
    show_default=True,
)
def long(
    ctx,
    genome,
    longreads,
    output,
    threads,
    database,
    force,
    bin_size,
    cutoff,
    width,
    min_orfs,
    phrog_fract,
    median_z_cutoff,
    sub_sample,
    coverage,
    no_extra,
    memory_limit,
    all_zscores,
    min_contig_len,
    join_window,
    expect_mad_zero,
    evalue,
    min_reads,
    **kwargs,
):
    """Runs hafeZ with ONT long reads"""

    # validates the directory  (need to before I start hafez or else no log file is written)
    instantiate_dirs(output, force)

    logdir: Path = Path(output) / "logs"

    output: Path = Path(output)

    params = {
        "--genome": genome,
        "--longreads": longreads,
        "--output": output,
        "--threads": threads,
        "--database": database,
        "--force": force,
        "--bin_size": bin_size,
        "--cutoff": cutoff,
        "--width": width,
        "--min_orfs": min_orfs,
        "--phrog_fract": phrog_fract,
        "--median_z_cutoff": median_z_cutoff,
        "--sub_sample": sub_sample,
        "--coverage": coverage,
        "--no_extra": no_extra,
        "--memory_limit": memory_limit,
        "--all_zscores": all_zscores,
        "--min_contig_len": min_contig_len,
        "--expect_mad_zero": expect_mad_zero,
        "--evalue": evalue,
        "--min_reads": min_reads,
    }

    # initial logging and list all the params
    start_time = begin_hafeZ(params)

    # validates fasta
    validate_fasta(genome)

    # check dependencies
    check_dependencies()

    # checks database installation
    logger.info(f"Checking database installation at {database}.")
    check_db_installation(database, install_flag=False)

    logger.info(f"Checking memory limit of {memory_limit}.")
    check_memory_limit(memory_limit)

    # instanatiate the premap class with coverage too
    ## premap.py
    premap = Premap()
    premap.coverage = coverage
    premap.output = output
    premap.start_time = start_time

    # filter the input genome to get long enough contigs
    premap.process_fasta(genome, min_contig_len)

    # get genome length of contigs
    premap.get_genome_length()

    # estimate coverage
    logger.info("Estimating overall genome coverage.")
    premap.calculate_estimated_coverage_long(longreads)
    logger.info(f"Estimated genome coverage is {premap.estimated_coverage:.2f}x.")

    # subsample if coverage is too high and set the Paths to the reads
    longreads = premap.subsample_long(longreads, logdir)

    #######################
    # mapping.py
    # no class needed
    ########################
    # run the mapping commands - no class as files written to temp
    logger.info("Mapping reads to filtered genome and sorting.")

    # map to the filtered FASTA
    seq_temp_out: Path = output / "temp_genome.fasta"
    minimap_long(longreads, seq_temp_out, output, threads, logdir)
    get_bam(output, threads, memory_limit, logdir)
    logger.info("Calculating per base coverage.")
    get_cov(output, threads, logdir)

    #######################
    # mapping_calcs.py
    #######################
    # smooths depths
    ## note it fails if the induction is too clean!
    # i.e. if you just simulate the phage reads, so need a little bit of background

    logger.info("Smoothing signals.")
    depths, cov = smooth_depths(output, bin_size, threads)

    # get Z scores
    zscores, median, mad = get_ZScores(depths, output, start_time, cov, expect_mad_zero)

    ########################
    # ROI
    # process_roi.py
    # uses Haf class to store everything
    ########################

    # instanatiate the Haf class
    haf = Haf()
    haf.all_zscores = all_zscores
    haf.median = median
    haf.mad = mad
    haf.zscores = zscores
    haf.output = output
    haf.depths = depths
    haf.start_time = start_time

    logger.info("Finding Regions of Interest (ROIs).")
    haf.get_ROIs(cutoff, width)
    logger.info("Joining any ROIs that are close to each other.")
    # RoIs within this window will be joined
    haf.join_ROIs(join_window)

    # filter rois based on width
    logger.info("Removing any small ROIs.")
    haf.filter_ROIs(width)

    ### differs for long from here

    logger.info("Extracting depth information for each ROI.")
    # number of bases before and after the ROI
    neighbourhood = 15001
    logger.info("Parsing and processing the BAM file.")
    haf.bamfile_to_pandas_long(output, neighbourhood)

    # find ROIs that span an entire contig
    logger.info("Checking for ROIs near contig ends.")
    haf.find_contig_end_rois(neighbourhood)

    #### LONG: find evidence of  reads that map twice - primary and supp - this means they must be over 2 different parts of the genome (start and end of prophage) ####

    # minimum number of multimapped reads that need to be on each side of the ROI for the phage to be induced
    min_reads = 5

    # find evidence of pairs of reads that are distant from each other (i.e. look to span deletion)
    logger.info("Finding reads that span both ends of the ROI.")
    haf.find_multimapped_reads_long(neighbourhood, min_reads)

    # adds relevant columns for downstream postprocessing

    roi_df = haf.roi_df
    roi_df["roi"] = roi_df["accession"]
    roi_df["start_pos"] = roi_df["likely_start"]
    roi_df["end_pos"] = roi_df["likely_end"]
    roi_df["circular"] = False
    roi_df["contig_split"] = "NaN"

    ###########################
    # orfs.py
    ###########################

    ##### HMMs
    # roi sequences
    logger.info("Predicting CDS on all ROIs using Pyrodigal.")
    roi_dna = get_roi_sequences(roi_df, premap.filtered_seq_dict, output)

    #### call orfs in genome ####
    # orf df, AA and nt
    orf_df, orfs_aa, orfs_dna = get_orfs(premap.filtered_seq_dict, premap.multicontig)

    ### extract orfs for roi ###

    roi_df, roi_orf_aa, roi_orf_dna = extract_roi_orfs(
        orf_df,
        orfs_aa,
        roi_df,
        orfs_dna,
        output,
        min_orfs,
        all_zscores,
        depths,
        median,
        mad,
        start_time,
    )

    #### Screen all roi orfs vs phrogs db with pyhmmer ####

    logger.info("Running PyHMMER on PHROGs.")
    best_phrog_results_dict = run_pyhmmer(database, output, threads, evalue=evalue)

    # calculate the percentage of orfs with PHROG hit for each roi
    logger.info("Calculating the proportion of ROI CDS with PHROGs hits.")
    roi_df = calc_phrogs_frac(
        best_phrog_results_dict,
        roi_df,
        phrog_fract,
        evalue,
        output,
        all_zscores,
        depths,
        median,
        mad,
        start_time,
    )

    #######################
    #### post_processing.py
    #######################

    # give rois new, final, names ahead of processing outputs

    roi_df = get_names(roi_df)

    # maybe can skip for long?
    logger.info("Finding possible attatchment sites.")
    roi_df = get_att(roi_df, premap.filtered_seq_dict, output, logdir)

    # save multifasta of roi genome seqs to file
    logger.info("Saving ROI sequences.")
    output_roi_seqs(roi_df, roi_dna, output)

    # output graph showing positions
    logger.info("Output plots.")
    output_prophage_graphs(roi_df, depths, output, median, mad)

    # output hmm table
    logger.info("Writing PyHMMER Output to file.")
    output_all_phrogs(roi_orf_aa, best_phrog_results_dict, output, database, evalue)

    #### output summary table of rois found ####
    logger.info("Writing ROI summary table.")
    output_roi_table(roi_df, output)

    # output roi orf aa and dna sequences
    logger.info("Writing ROI AA and DNA sequences.")
    output_roi_orfs(roi_orf_dna, roi_df, output, roi_orf_aa)

    # make Z-score graphs for each contig if -Z flag given
    if all_zscores is True:
        logger.info("Making Z-score graphs.")
        output_contig_Z(depths, output, median, mad)

    # exit hafeZ
    exit_success(output, start_time)


"""
database command
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@click.option(
    "-d",
    "--database",
    help="Specific path to download and install the hafeZ database directory. Optional. ",
    default=None,
    show_default=True,
)
def database(ctx, database):
    """Downloads and installs hafeZ database"""
    logger.add(lambda _: sys.exit(1), level="ERROR")
    # set defaultdefault
    if database is None:
        current_dir = os.path.dirname(os.path.realpath(__file__))
        database = os.path.join(current_dir, "../", "databases/")
        # make the database includes recursive
        if not os.path.exists(database):
            os.makedirs(database)

    # download the db
    check_db_installation(database, install_flag=True)


@click.command()
def citation(**kwargs):
    """Print the citation for hafeZ"""
    print_citation()


main_cli.add_command(short)
main_cli.add_command(long)
main_cli.add_command(database)
main_cli.add_command(citation)


def main():
    main_cli()


if __name__ == "__main__":
    main()
