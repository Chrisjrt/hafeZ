#!/usr/bin/env python3

##############################################
########### IMPORT ALL LIBRARIES #############
##############################################
import argparse
import os
import sys
import time
import shutil
import hZ.get_db
import hZ.process_genome
import hZ.mapping
import hZ.mapping_calcs
import hZ.process_rois
import hZ.get_output
import pandas as pd
import numpy as np
##############################################



##############################################
############ SET CURRENT VERSION #############
##############################################

__version__ = '1.0.3'

##############################################



##############################################
################ BEGIN MAIN ##################
##############################################

def main():
    '''Run hafeZ using command line arguments'''

    ##############################################
    ## CREATE ARGUMENT PARSER AND SET ARGUMENTS ##
    ##############################################

    #### create parser ####
    parser = argparse.ArgumentParser(
                prog='hafeZ',
                description='Identify inducible prophages through bacterial genomic read mapping. Minimum required input outlined above.',
                usage = ' setup:  hafeZ.py -G db_path\n\trun: hafeZ.py -f assembly.fasta -r1 reads.fastq.gz -r2 reads.fastq.gz -o output_folder -D db_folder',
                add_help=False)


    #### add setup arguments ####
    setup = parser.add_argument_group('required arguments for initial setup')
    # option to get/format pVOGs DB and save it to the given directory
    setup.add_argument(
        '-G', '--get_db',
        metavar = 'path',
        nargs = '?',
        help = 'use this option to get and format pVOGs database in the given directory',
        required = False)
    setup.add_argument(
        '-T', '--db_type',
        metavar = 'db',
        nargs = '?',
        help = 'choose which database you want to download, currently available ones are pVOGs or PHROGs',
        required = False)

    #### add arguments required for running hafeZ ####
    required = parser.add_argument_group('required arguments for running hafeZ')

    # assembly fasta input option
    required.add_argument(
        '-f', '--fasta',
        metavar = 'path',
        help = 'path to genome assembly in fasta format',
        nargs = '?',
        required =False)
    # read 1 input option
    required.add_argument(
         '-r1','--reads1',
         metavar = 'path',
         help = 'path to read set in fastq/fastq.gz format.',
         nargs = '?',
         required =False)
    # read 2 input option
    required.add_argument(
         '-r2','--reads2',
         metavar = 'path',
         help = 'path to second read set in fastq/fastq.gz format',
         nargs = '?',
         required =False)
    # output path input option
    required.add_argument(
         '-o','--output_folder',
         metavar = 'path',
         help = 'desired output folder path',
         nargs = '?',
         required =False)
    # option to give path to folder containing formatted pVOGs DB files if theyve already been downloaded
    required.add_argument(
        '-D', '--db_path',
        metavar = 'path',
        nargs = '?',
        help = "path to the directory containing the pVOGs files",
        required = False)

    #### add optional arguments ####
    optional = parser.add_argument_group('optional arguments')

    # bin size input option
    optional.add_argument(
        '-b', '--bin_size',
        metavar = 'int',
        help = 'set bin size in bp to use for coverage depth smoothing. Must be an odd number. (default = 3001)',
        nargs = '?',
        type = int,
        default = 3001,
        required =False)
    # cutoff input option
    optional.add_argument(
        '-c', '--cutoff',
        metavar = 'float',
        help = 'set Z-score cutoff for initially detecting rois (default = 3.5)',
        nargs = '?',
        type=float,
        default = 3.5,
        required =False)
    # width of roi input option
    optional.add_argument(
        '-w', '--width',
        metavar = 'int',
        help = 'set minimum width (bp) of roi that passes Z-score cutoff for initially detecting rois. (default = 4000)',
        nargs = '?',
        type = int,
        default = 4000,
        required =False)
    # min number of orfs in roi option
    optional.add_argument(
        '-m', '--min_orfs',
        metavar = 'int',
        help = 'set minimum number of ORFs needed in an roi (default = 6)',
        nargs = '?',
        type = int,
        default = 6,
        required =False)
    # min fraction of orfs hit by pvog option
    optional.add_argument(
        '-p', '--pvog_fract',
        metavar = 'float',
        help = 'set minimum fraction of ORFs in an roi showing homology to pVOGs (default = 0.1)',
        nargs = '?',
        type=float,
        default = 0.1,
        required =False)
    # min median Z-score of an roi to be kept option
    optional.add_argument(
        '-z', '--median_z_cutoff',
        metavar = 'float',
        help = 'set minimum median Z-score for an roi to be retained (default = 3.5)',
        nargs = '?',
        type=float,
        default = 3.5,
        required =False)
    # option to make zscore graphs of all contigs
    optional.add_argument(
        '-k', '--keep_threshold',
        metavar = 'int',
        help = 'set threshold for number of best soft clip combinations to keep for each roi (default=50) ',
        nargs = '?',
        type=float,
        default = 50,
        required =False
        )
    optional.add_argument(
        '-S', '--sub_sample',
        action='store_true',
        help = 'Randomly sub-sample reads to adjust overall mapping coverage of genome. N.B. Use -C/--coverage to adjust coverage desired',
        required =False)
    # desired coverage of genome to be used for subsampling reads
    optional.add_argument(
        '-C', '--coverage',
        metavar = 'float',
        help = 'set desired coverage of genome to be used for subsampling reads (default = 100.0). N.B. Must be used with -S\--sub_sample flag',
        nargs = '?',
        type=float,
        default = 100.0,
        required =False
        )
    # add argument to skip accuracy checks using clipped sequences
    optional.add_argument(
        '-N', '--no_extra',
        action='store_true',
        help = 'use to turn off extra accuracy checks using clipped sequences, might give same results, might give extra rois',
        required =False)
    # threads input option
    optional.add_argument(
         '-t', '--threads',
         metavar = 'int',
         help = 'set number of threads to use (default = 1)',
         nargs = '?',
         type = int,
         default = 1,
         required =False)
    optional.add_argument(
        '-M', '--memory_limit',
        metavar = 'limit',
        nargs = '?',
        default = '768M',
        help = "set upper bound per thread memory limit for samtools, suffix K/M/G recognized (default = 768M)",
        required = False)
    # option to allow output_folder to be overwritten
    optional.add_argument(
        '-O', '--overwrite',
        action='store_true',
        help = 'force overwrite of ouput folder if already exists',
        required =False)
    # option to make zscore graphs of all contigs
    optional.add_argument(
        '-Z', '--all_Zscores',
        action='store_true',
        help = 'make graphs of all contig Z-scores even if no roi found (useful for manual inspection). N.B. This will make graphs for each contig, so if you have 100 contigs you will get 100 graphs.',
        required =False)
    # get help option
    optional.add_argument('-h', '--help',
        action = 'help',
        help='show this help message and exit')
    # get version option
    optional.add_argument(
         '-v','--version',
         action='version',
         version = '%(prog)s '+__version__)


    ##### Activate arguments #####
    args = parser.parse_args()

    print('\nRunning hafeZ version {} with the following settings:\nhafeZ.py -f {} -r1 {} -r2 {} -o {} -D {} -c {} -b {} -w {} -m {} -t {} -p {} -z {}'.format(
        __version__, args.fasta, args.reads1, args.reads2, args.output_folder, args.db_path, args.cutoff, args.bin_size, args.width, args.min_orfs, args.threads, args.pvog_fract, args.median_z_cutoff))


    ##############################################
    ############ CATCH ARGUMENT ERRORS ###########
    ##############################################

    # check that dbs have been set properly:
    if args.get_db is not None and args.db_path is not None:
        parser.print_help()
        print('hafeZ: error: --get_db and --db_path cannot be used together, use one or the other')
        sys.exit(0)
    if args.get_db is None and args.db_path is None:
        parser.print_help()
        print('hafeZ: error: either --get_db or --db_path options must be used and should contain a path')
        sys.exit(0)

    if args.db_type is None or args.db_type.lower() not in ['pvogs','phrogs']:
        parser.print_help()
        print('hafeZ: error: desired --db_type must be given, currently it hasnt been set properly')
        sys.exit(0)

     # check memory limit has been given properly
    args.memory_limit = args.memory_limit.upper()
    if not any(x in args.memory_limit for x in ['K','G','M']):
        parser.print_help()
        print('hafeZ: error: memory limit format must be in format of e.g. 2GB or 2MB (GB being gigabytes of RAM and MB being megabytes of RAM)')
        sys.exit(1)


    ##############################################


    ##############################################
    ######## BEGIN RUNNING MAIN COMMANDS #########
    ##############################################

    print('\n{:#^50}'.format(''))
    print('{:#^50}'.format(' Running hafeZ '))
    print('{:#^50}'.format(''))

    run_start_time = time.time()

    #### get pvogs db if needed ####

    if args.get_db != None:
        if not os.path.exists(args.get_db):
            os.makedirs(args.get_db, exist_ok=False)
        if args.db_type.lower() == "pvogs":
            hZ.get_db.get_dbs_pvogs(args.get_db)
            hZ.get_db.process_vogsTable_pvogs(args.get_db)
            hZ.get_db.process_vogsHMM_pvogs(args.get_db)
        elif args.db_type.lower() == "phrogs":
            hZ.get_db.get_dbs_phrogs(args.get_db)
            # hZ.get_db.process_vogsTable_phrogs(args.get_db)
            # hZ.get_db.process_vogsHMM_phrogs(args.get_db)
        else:
            parser.print_help()
            print('hafeZ: error: --db_type given is not recognised, see help for available databases')
            sys.exit(0)

        args.db_path = args.get_db

    #### check that all arguments have been given for running

    if (args.output_folder != None and
        args.fasta != None and
        args.reads1 != None and
   	args.reads2 != None):

        #### make results output directory ####

        if args.overwrite == True:
            if os.path.exists(args.output_folder):
                shutil.rmtree(args.output_folder)
            os.makedirs(args.output_folder, exist_ok=True)
        else:
            os.makedirs(args.output_folder, exist_ok=True)
        args.output_folder = os.path.abspath(args.output_folder)

        #### check to see if multicontig + filter fasta to remove small contigs ####


        seq, multicontig, fasta_filtered, raw_seq_dict = hZ.process_genome.process_fasta(args.fasta,args.output_folder)
        if isinstance(seq, str):
            hZ.get_output.output_no_roi(args.output_folder)
            hZ.get_output.clean_up(args.output_folder)
            run_end_time = '{:.2f}'.format(time.time() - run_start_time)
            print('\n{:#^50}'.format(''))
            print('{:#^50}'.format(' hafeZ run complete! '))
            print('{:#^50}'.format(' Total run time: ' + run_end_time + ' seconds '))
            print('{:#^50}'.format(''))
            sys.exit(0)

        #### get overall length of genome ####

        genome_length = hZ.process_genome.get_genome_length(raw_seq_dict)


        if args.sub_sample == True:
            #### caluclate rough estimated coverage ####

            coverage = hZ.mapping.get_rough_coverage_prediction(args.reads1, args.reads2, genome_length)


            if coverage > args.coverage:
                #### take random subsample of reads to appropriate coverage level

                args.reads1, args.reads2 = hZ.mapping.take_read_sample(args.reads1, args.reads2, coverage, args.output_folder, args.coverage)
            else:
                print('True coverage ({}x) is less than desired coverage ({}x).\nContinuing with True coverage.'.format(round(coverage,1), round(args.coverage,1)))


        #### map reads to assembly ####

        hZ.mapping.minimap_paired(args.reads1, args.reads2, fasta_filtered, args.output_folder, args.threads)


        #### convert sam to bam ####

        hZ.mapping.get_bam(args.output_folder, args.threads, args.memory_limit)


        #### index and get coverage depths ####

        hZ.mapping.get_cov(args.output_folder,args.threads)


        #### smooth the depths ####

        depths, cov = hZ.mapping_calcs.smooth_depths(args.output_folder, args.bin_size, args.threads)


        #### get Z scores ####

        zscores, median, mad = hZ.mapping_calcs.get_ZScores(depths, multicontig)


        #### get roi coords and begin building output df ####

        roi_df = hZ.mapping_calcs.get_ROIs(zscores,args.cutoff,args.output_folder,args.width)
        if isinstance(roi_df, str):
            hZ.get_output.output_no_roi(args.output_folder)
            hZ.get_output.clean_up(args.output_folder)
            if args.all_Zscores == True:
                hZ.get_output.output_contig_Z(depths,args.output_folder,median, mad)
            run_end_time = '{:.2f}'.format(time.time() - run_start_time)
            print('\n{:#^50}'.format(''))
            print('{:#^50}'.format(' hafeZ run complete! '))
            print('{:#^50}'.format(' Total run time: ' + run_end_time + ' seconds '))
            print('{:#^50}'.format(''))
            sys.exit(0)


        #### join rois that are near each other ####

        roi_df = hZ.process_rois.join_ROIs(roi_df)


        #### filter rois based on width ####

        roi_df = hZ.process_rois.filter_ROIs(roi_df,args.width)
        if isinstance(roi_df, str):
            hZ.get_output.output_no_roi(args.output_folder)
            hZ.get_output.clean_up(args.output_folder)
            if args.all_Zscores == True:
                hZ.get_output.output_contig_Z(depths,args.output_folder,median, mad)
            run_end_time = '{:.2f}'.format(time.time() - run_start_time)
            print('\n{:#^50}'.format(''))
            print('{:#^50}'.format(' hafeZ run complete! '))
            print('{:#^50}'.format(' Total run time: ' + run_end_time + ' seconds '))
            print('{:#^50}'.format(''))
            sys.exit(0)


        #### process sam file ####

        sam_df = hZ.process_rois.samfile_to_pandas(args.output_folder, roi_df, args.threads)


        #### find any rois that span an entire contig ####

        end_roi_df, roi_df = hZ.process_rois.find_contig_end_rois(roi_df)

        #### find evidence of pairs of reads that are distant from each other (i.e. look to span deletion) ####

        if not isinstance(roi_df, str):
            roi_df = hZ.process_rois.find_far_reads(sam_df, args.width, roi_df) ##### might want to make this a non exclusion stage and instead an informative stage (i.e. dont delete those that dont have, instead just make column value that says it )

        #### combine roi_df and end_df for next steps: ####
        if not isinstance(roi_df, str) and not isinstance(end_roi_df, str):
            roi_df = pd.concat([roi_df, end_roi_df])
            roi_df = roi_df.reset_index(drop=True)
            end_roi_df = 'exit'
        elif isinstance(roi_df, str) and not isinstance(end_roi_df, str):
            roi_df = end_roi_df
            end_roi_df = 'exit'

        #### find soft clipped reads  ####

        if not isinstance(roi_df, str):
            roi_df, clips_df = hZ.process_rois.find_soft_clippings(sam_df,roi_df,args.width, args.threads)#, fasta_filtered,args.output_folder, args.threads)


        #### refind rois at ends ####

        if not isinstance(roi_df, str):
            roi_df = hZ.process_rois.find_contig_end_rois_again(roi_df)


        #### calculate median Z-scores and filter rois ####

        if not isinstance(roi_df,str):
            #### calculate median Z-score for each roi ####
            roi_df = hZ.process_rois.calc_roi_Z_medians(roi_df,cov,median,mad,args.median_z_cutoff, args.threads)
            #### filter rois using median Z-score for each roi ####
            roi_df = hZ.process_rois.filter_Z_medians(roi_df,args.median_z_cutoff)
            ### filter to keep only best X rois ####
            if not isinstance(roi_df,str):
                roi_df = hZ.process_rois.keep_only_x_best(roi_df,args.keep_threshold)

        #### check to see if any rois are circular ####
        if not isinstance(roi_df, str):
            circular_df, roi_df = hZ.process_rois.check_for_circular_roi(sam_df,roi_df)#,args.threads)
        else:
            circular_df = 'exit'


        #### do the extra accuracy steps ####

        if args.no_extra == False:
            ### collect the ends of the soft clipped reads found, map them, collect the results, and process them ####
            if not isinstance(roi_df, str):
                roi_df = hZ.process_rois.collecting_clipped_reads(roi_df, clips_df, args.output_folder,args.threads)
                hZ.mapping.map_clipped_reads(args.output_folder,args.threads,args.fasta,  args.memory_limit)
                roi_df, end_roi_df = hZ.process_rois.seperate_rois_near_ends(roi_df)
            if not isinstance(roi_df, str):
                clip_df = hZ.process_rois.process_clip_sam(args.output_folder,roi_df)
                if clip_df.empty:
                    roi_df = 'exit'
                else:
                    roi_df = hZ.process_rois.get_clip_pos(clip_df, roi_df, args.output_folder)
            if not isinstance(end_roi_df, str):
                clip_df2 = hZ.process_rois.process_clip_sam_end_rois(args.output_folder, end_roi_df)
                end_roi_df = hZ.process_rois.get_clip_pos_end_rois(clip_df2, end_roi_df, raw_seq_dict)
            if not isinstance(end_roi_df, str):
                end_roi_df = hZ.process_rois.filter_by_end_status(end_roi_df)
            if not isinstance(end_roi_df, str):
                end_roi_df = hZ.process_rois.reformat_end_roi_tables(end_roi_df, raw_seq_dict)


        #### rejoin all databases together ####
        df_list = []
        for i in [roi_df,end_roi_df,circular_df]:
            if not isinstance(i, str):
                df_list.append(i)
        if len(df_list) > 1:
            roi_df = pd.concat(df_list)
            roi_df = roi_df.reset_index(drop=True)
        elif len(df_list) == 1:
            roi_df = df_list[0]
        elif len(df_list) == 0:
            hZ.get_output.output_no_roi(args.output_folder)
            hZ.get_output.clean_up(args.output_folder)
            if args.all_Zscores == True:
                hZ.get_output.output_contig_Z(depths,args.output_folder,median, mad)
            run_end_time = '{:.2f}'.format(time.time() - run_start_time)
            print('\n{:#^50}'.format(''))
            print('{:#^50}'.format(' hafeZ run complete! '))
            print('{:#^50}'.format(' Total run time: ' + run_end_time + ' seconds '))
            print('{:#^50}'.format(''))
            sys.exit(0)

        # pick best rois if read checking has been done
        if args.no_extra == True:
            roi_df = hZ.process_rois.quick_filter(roi_df)

        #### extract roi sequences ####
        for index, row in roi_df.iterrows():  #### move this elswhere
            if pd.isna(row['circular']):
                roi_df.loc[index, 'circular'] = False
            else:
                roi_df.loc[index,'med_z'] = np.nan

        ### this may be a good spot to filter the circular rois by z score

        roi_dna = hZ.process_genome.get_roi_sequences(roi_df, seq, args.output_folder)

        #### call orfs in genome ####

        orf_df, orfs_aa, orfs_dna = hZ.process_genome.get_orfs(seq, multicontig)


        ### extract orfs for roi ###

        roi_df, roi_orf_aa, roi_orf_dna = hZ.process_rois.extract_roi_orfs(orf_df, orfs_aa, roi_df, orfs_dna, args.output_folder, args.min_orfs)
        if isinstance(roi_df, str):
            hZ.get_output.output_no_roi(args.output_folder)
            hZ.get_output.clean_up(args.output_folder)
            if args.all_Zscores == True:
                hZ.get_output.output_contig_Z(depths,args.output_folder,median, mad)
            run_end_time = '{:.2f}'.format(time.time() - run_start_time)
            print('\n{:#^50}'.format(''))
            print('{:#^50}'.format(' hafeZ run complete! '))
            print('{:#^50}'.format(' Total run time: ' + run_end_time + ' seconds '))
            print('{:#^50}'.format(''))
            sys.exit(0)

        if args.db_type.lower() == 'pvogs':
            ### Screen all roi genes vs pvogs db ####
            hmm_df = hZ.process_rois.hmm_scan(args.output_folder, args.threads, args.db_path)
            ### calculate fraction of orfs hit by pVOGs ####
            roi_df = hZ.process_rois.calc_pVOGs_frac(hmm_df,roi_df,args.pvog_fract)
            if isinstance(roi_df, str):
                hZ.get_output.output_no_roi(args.output_folder)
                hZ.get_output.clean_up(args.output_folder)
                if args.all_Zscores == True:
                    hZ.get_output.output_contig_Z(depths,args.output_folder,median, mad)
                run_end_time = '{:.2f}'.format(time.time() - run_start_time)
                print('\n{:#^50}'.format(''))
                print('{:#^50}'.format(' hafeZ run complete! '))
                print('{:#^50}'.format(' Total run time: ' + run_end_time + ' seconds '))
                print('{:#^50}'.format(''))
                sys.exit(0)
        elif args.db_type.lower() == 'phrogs':
            #### Screen all roi genes vs phrogs db ####
            hmm_df = hZ.process_rois.hhblits_phrogs(args.output_folder, args.threads, args.db_path)
            #### calculate fraction of orfs hit by phrogs ####
            roi_df = hZ.process_rois.calc_phrogs_frac(hmm_df,roi_df,args.pvog_fract)
            if isinstance(roi_df, str):
                hZ.get_output.output_no_roi(args.output_folder)
                hZ.get_output.clean_up(args.output_folder)
                if args.all_Zscores == True:
                    hZ.get_output.output_contig_Z(depths,args.output_folder,median, mad)
                run_end_time = '{:.2f}'.format(time.time() - run_start_time)
                print('\n{:#^50}'.format(''))
                print('{:#^50}'.format(' hafeZ run complete! '))
                print('{:#^50}'.format(' Total run time: ' + run_end_time + ' seconds '))
                print('{:#^50}'.format(''))
                sys.exit(0)

        ### give rois new, final, names ahead of processing outputs ####

        roi_df = hZ.get_output.get_names(roi_df)

        ### find atts ####
        roi_df = hZ.process_rois.get_att(roi_df,seq,args.output_folder)

        ### output multifasta of roi genome seqs ####

        hZ.get_output.output_roi_seqs(roi_df,roi_dna,args.output_folder)

        #### output graph showing positions ####

        hZ.get_output.output_prophage_graphs(roi_df,depths,args.output_folder, median, mad)


        #### fix hmm table ####

        if args.db_type.lower() == 'pvogs':
            hZ.get_output.output_hmm_table(roi_df,args.output_folder)
        elif args.db_type.lower() == 'phrogs':
            hZ.get_output.output_hmm_table_phrogs(roi_df,args.output_folder,args.db_path)

        #### output summary table of rois found ####

        hZ.get_output.output_roi_table(roi_df,args.output_folder)

        #### output roi orf aa and dna sequences ####

        hZ.get_output.output_roi_orfs(roi_orf_dna,roi_df,args.output_folder,roi_orf_aa)

        #### make Z-score graphs for each contig if -Z flag given ####

        if args.all_Zscores == True:
            hZ.get_output.output_contig_Z(depths,args.output_folder,median, mad)


        #### tidy up by removing temp files ####

        hZ.get_output.clean_up(args.output_folder)


        #### print end time ####

        run_end_time = '{:.2f}'.format(time.time() - run_start_time)

        print('\n{:#^50}'.format(''))
        print('{:#^50}'.format(' hafeZ run complete! '))
        print('{:#^50}'.format(' Total run time: ' + run_end_time + ' seconds '))
        print('{:#^50}'.format(''))

######################################

if __name__ == '__main__' :
    main()
