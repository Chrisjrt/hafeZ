import time
from Bio import SeqIO
import numpy as np
import subprocess
import pandas as pd
from io import StringIO
import matplotlib.pyplot as plt
import sys
import re
from scipy import stats
import pysam
from collections import Counter
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import multiprocessing
from functools import partial
import hZ.misc_functions
import itertools
import os

def join_ROIs(roi_df):
    print('\n{:#^50}'.format(' Merging close regions of interest '))
    start_time = time.time()
    counter_1 = 0
    joined_roi_df = {'accession': [], 'start': [], 'end': [], 'length': [], 'contig_len': [], 'contig': []}
    for index, row in roi_df.iterrows():
        counter_1 = counter_1 + 1
        past_df_len = 0
        sub_df = roi_df[(roi_df['contig'] == row['contig']) &
                               (((roi_df['start'].between(row['start'] - 10000,row['end'] + 10000)) | (roi_df['end'].between(row['start'] - 10000,row['end'] + 10000))) |
                               (((row['start'] >= roi_df['start'] - 10000) & (row['start']<= roi_df['end'] + 10000)) | ((row['end'] >= roi_df['start'] - 10000) & (row['end'] <= roi_df['end'] + 10000))))]
        current_df_len = len(sub_df)
        sub_list = list(sub_df['start']) + list(sub_df['end'])
        min_pos = min(sub_list)
        max_pos = max(sub_list)
        while current_df_len > past_df_len:
            past_df_len = current_df_len
            sub_df = roi_df[(roi_df['contig'] == row['contig']) &
                                ((roi_df['start'].between(min_pos - 10000, max_pos + 10000) | roi_df['end'].between(min_pos - 10000, max_pos + 10000)) |
                                (((min_pos >= roi_df['start'] - 10000) & (min_pos <= roi_df['end'] + 10000)) | ((max_pos >= roi_df['start'] - 10000) & (max_pos <= roi_df['end'] + 10000))))]
            sub_list = list(sub_df['start']) + list(sub_df['end'])
            min_pos = min(sub_list)
            max_pos = max(sub_list)
            current_df_len = len(sub_df)
        joined_roi_df['accession'].append('~'.join(sub_df['accession'].iloc[0].split('~')[:-1]) + '~' + str(counter_1) )
        joined_roi_df['start'].append(min_pos)
        joined_roi_df['end'].append(max_pos)
        joined_roi_df['length'].append(max_pos - min_pos)
        joined_roi_df['contig_len'].append(sub_df['contig_len'].iloc[0])
        joined_roi_df['contig'].append(sub_df['contig'].iloc[0])
    joined_roi_df = pd.DataFrame.from_dict(joined_roi_df)
    joined_roi_df = joined_roi_df.drop_duplicates(subset=['start','end','length','contig_len','contig'],keep='first')
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return joined_roi_df

def filter_ROIs(roi_df,width):
    print('\n{:#^50}'.format(' Removing small regions of interest '))
    start_time = time.time()
    roi_df = roi_df[roi_df['length'] >= width]
    if len(roi_df) < 1:
        end_time = '{:.4f}'.format(time.time() - start_time)
        print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
        return 'exit'
    end_time = '{:.4f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return roi_df

def samfile_to_pandas(output_folder, roi_df, threads):
    print('\n{:#^50}'.format(' Parsing and processing sam file '))
    start_time = time.time()
    samfile = pysam.AlignmentFile(output_folder + '/temp_minimap_sorted.bam', 'rb')
    df_dict = {'qname':[], 'rname':[], 'pos':[], 'cigar':[],'tlen':[],'seq':[],'flag':[],'pnext':[], 'qqual':[]}
    for index, row in roi_df.iterrows():
        contig = row['accession'].split('~')[0]
        start = row['start'] - 15001
        if start < 0:
            start = 0
        end = row['end'] + 15001
        if end > row['contig_len']:
            end = row['contig_len']
        for reads in samfile.fetch(contig, start, end):
            df_dict['qname'].append(reads.query_name)
            df_dict['rname'].append(contig)
            df_dict['pos'].append(reads.reference_start)
            df_dict['cigar'].append(reads.cigarstring)
            df_dict['tlen'].append(reads.template_length)
            df_dict['seq'].append(reads.query_sequence)
            df_dict['flag'].append(reads.flag)
            df_dict['pnext'].append(reads.next_reference_start)
            df_dict['qqual'].append(reads.query_qualities.tolist())
    df = pd.DataFrame.from_dict(df_dict)
    df = df.drop_duplicates(subset = ['qname', 'rname', 'pos', 'cigar', 'tlen', 'seq', 'flag', 'pnext'], keep='first')
    df['tlen'] = df['tlen'].astype(int)
    df['pos'] = df['pos'].astype(int)
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return df

def find_contig_end_rois(roi_df):
    print('\n{:#^50}'.format(' Checking for rois near contig ends '))
    start_time = time.time()
    for index, row in roi_df.iterrows():
        roi_df.loc[index, 'center'] = (row['start'] + row['end'])/2
    end_df = roi_df[((roi_df['start'].between(0,15001)) | (roi_df['end'].between(roi_df['contig_len'] - 15001, roi_df['contig_len'])))].copy() # pull out rois into a seperate df that would have failed later tests unfairly due to being close to end of contig
    roi_df = roi_df[~((roi_df['start'].between(0,15001)) | (roi_df['end'].between(roi_df['contig_len'] - 15001, roi_df['contig_len'])))].copy()
    for index, row in end_df.iterrows():
        if (0 <= row['start'] <= 15001) and (row['contig_len'] - 15001 <= row['end'] <= row['contig_len']):
            end_df.loc[index, 'end_info'] = 'whole'
        elif (0 <= row['start'] <= 15001) and (not row['contig_len'] - 15001 <= row['end'] <= row['contig_len']):
            end_df.loc[index, 'end_info'] = 'left'
        elif (not 0 <= row['start'] <= 15001) and (row['contig_len'] - 15001 <= row['end'] <= row['contig_len']):
            end_df.loc[index, 'end_info'] = 'right'
    roi_df['end_info'] = np.nan
    if len(roi_df) < 1:
        roi_df = 'exit'
    if len(end_df) < 1:
        end_df = 'exit'
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return end_df, roi_df

def find_contig_end_rois_again(roi_df):
    print('\n{:#^50}'.format(' Checking for rois near contig ends again '))
    start_time = time.time()
    for index, row in roi_df.iterrows():
        roi_df.loc[index, 'center'] = (row['start_pos'] + row['end_pos'])/2
    end_df = roi_df[((roi_df['start_pos'].between(0,15001)) | (roi_df['end_pos'].between(roi_df['contig_len'] - 15001, roi_df['contig_len'])))].copy() # pull out rois into a seperate df that would have failed later tests unfairly due to being close to end of contig
    roi_df = roi_df[~((roi_df['start_pos'].between(0,15001)) | (roi_df['end_pos'].between(roi_df['contig_len'] - 15001, roi_df['contig_len'])))].copy()
    for index, row in end_df.iterrows():
        if (0 <= row['start_pos'] <= 15001) and (row['contig_len'] - 15001 <= row['end_pos'] <= row['contig_len']):
            end_df.loc[index, 'end_info'] = 'whole'
        elif (0 <= row['start_pos'] <= 15001) and (not row['contig_len'] - 15001 <= row['end_pos'] <= row['contig_len']):
            end_df.loc[index, 'end_info'] = 'left'
        elif (not 0 <= row['start_pos'] <= 15001) and (row['contig_len'] - 15001 <= row['end_pos'] <= row['contig_len']):
            end_df.loc[index, 'end_info'] = 'right'
    roi_df['end_info'] = np.nan
    roi_df = pd.concat([roi_df,end_df])
    roi_df = roi_df.reset_index(drop=True)
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return roi_df

def find_far_reads(sam_df, width, roi_df):
    print('\n{:#^50}'.format(' Finding distant read pairs '))
    start_time = time.time()
    df = sam_df[((sam_df['tlen'] > width) | (sam_df['tlen'] < (width*-1))) & (sam_df.duplicated(subset=['qname'], keep=False))].copy()
    for index, row in roi_df.iterrows():
        center = row['center']
        sub_df = df[(df['rname'] == row['contig']) &
                    (((df['pos'].between(row['start'] - 15000, center)) & (df['pnext'].between(center, row['end'] + 15000))) |
                    ((df['pnext'].between(row['start'] - 15000, center)) & (df['pos'].between(center, row['end'] + 15000))))].copy()
        sub_df = sub_df[sub_df.duplicated(subset=['qname'], keep=False)]
        if (len(sub_df)/2) < 1:
            roi_df = roi_df[roi_df['accession'] != row['accession']].copy()
        else:
            roi_df.loc[index,'far_start'] = np.median(sub_df['pos'][sub_df['pos'] < center].copy())
            roi_df.loc[index,'far_end'] = np.median(sub_df['pos'][sub_df['pos'] > center].copy())
            roi_df.loc[index,'far_count_pairs'] = len(sub_df.copy())/2
    if len(roi_df) < 1:
        end_time = '{:.2f}'.format(time.time() - start_time)
        print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
        return 'exit'
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return roi_df

def multiprocess_find_soft_clippings(roi_df, clips, width):
    clip_dict = {'roi':[], 'start_pos':[], 'start_count':[], 'end_pos':[], 'end_count':[], 'length':[], 'contig_len':[], 'contig':[], 'end_info':[]}
    for index, row in roi_df.iterrows():
        s_list = []
        sub_df = clips[(clips['rname'] == row['contig']) & (clips['pos'].between(row['start'] - 15000, row['end'] + 15000))]
        for i, r in sub_df.iterrows():
            M = re.sub("[0-9]+", "",r['cigar']).find('M')
            S = re.sub("[0-9]+", "",r['cigar']).find('S')
            if S == 1:
                s_list.append(r['pos'] + int(re.split(r'\D+',r['cigar'])[M]))
            else:
                s_list.append(r['pos'])
        if len(s_list) > 1:
            sub_df = pd.Series(s_list).value_counts().rename_axis('coords').reset_index(name='counts')
            sub_df = sub_df[sub_df['counts'] >= 10]
        if (len(sub_df[sub_df['coords'] < row['center']]) > 0) and (len(sub_df[sub_df['coords'] > row['center']]) > 0):
            for i, r in sub_df[sub_df['coords'] < row['center']].iterrows():
                for ind, rw in sub_df[sub_df['coords'] > row['center']].iterrows():
                    clip_dict['roi'].append(row['accession'])
                    clip_dict['start_pos'].append(r['coords'])
                    clip_dict['start_count'].append(r['counts'])
                    clip_dict['end_pos'].append(rw['coords'])
                    clip_dict['end_count'].append(rw['counts'])
                    clip_dict['length'].append(rw['coords'] - r['coords'])
                    clip_dict['contig_len'].append(row['contig_len'])
                    clip_dict['contig'].append(row['contig'])
                    clip_dict['end_info'].append(row['end_info'])
    df = pd.DataFrame.from_dict(clip_dict)
    df = df[df['length'] >= width]
    return df

def find_soft_clippings(sam_df,roi_df,width,threads):
    print('\n{:#^50}'.format(' Finding clipped reads '))
    start_time = time.time()
    clips = sam_df[(sam_df['cigar'].str.count('M') == 1) & (sam_df['cigar'].str.count('S') == 1) & (sam_df['cigar'].str.count('I') == 0) & (sam_df['cigar'].str.count('X') == 0) & (sam_df['cigar'].str.count('D') == 0)]
    chunks, pool = hZ.misc_functions.split_df_for_multiproc(roi_df,threads)
    df = pd.concat(pool.map(partial(multiprocess_find_soft_clippings, clips = clips, width=width), chunks))
    df = df.reset_index(drop=True)
    if len(df) < 1:
        return 'exit', 'exit'
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return df, clips

def collect_clips(roi_df, clips):
    read_list = []
    for index, row in roi_df.iterrows():
        sub_df = clips[(clips['rname'] == row['contig']) & (clips['pos'].between(row['start_pos'] - 200, row['end_pos'] + 200))]
        counter_2 = 1
        for i, r in sub_df.iterrows():
            M = re.sub("[0-9]+", "",r['cigar']).find('M')
            S = re.sub("[0-9]+", "",r['cigar']).find('S')
            if S == 1:
                matched = int(re.split(r'\D+',r['cigar'])[M])
                pos = r['pos'] + matched
                clipped_seq = r['seq'][matched:]
                clipped_qual = r['qqual'][matched:]
            else:
                matched = int(re.split(r'\D+',r['cigar'])[S])
                pos = r['pos']
                clipped_seq = r['seq'][:matched]
                clipped_qual = r['qqual'][:matched]
            if (pos == row['start_pos']) or (pos == row['end_pos']):
                clip_name = row['roi'] + '~' + str(pos) + '__' + str(counter_2)
                clipped = SeqRecord(Seq(clipped_seq), id=clip_name, name=clip_name, description='', dbxrefs=[])
                clipped.letter_annotations["phred_quality"] = clipped_qual
                read_list.append(clipped)
                counter_2 = counter_2 + 1
    return [read_list, roi_df]

def collecting_clipped_reads(roi_df, clips, output_folder, threads):
    print('\n{:#^50}'.format(' Collecting clipped ends '))
    start_time = time.time()
    roi_df = roi_df.reset_index(drop=True)
    chunks, pool = hZ.misc_functions.split_df_for_multiproc(roi_df,threads)
    values = pool.map(partial(collect_clips, clips = clips), chunks)
    read_list = list(itertools.chain.from_iterable([item[0] for item in values]))
    roi_df = pd.concat([item[1] for item in values])
    SeqIO.write(read_list, output_folder + '/temp_soft_clipped.fastq', 'fastq')
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return roi_df

def seperate_rois_near_ends(roi_df):
    print('\n{:#^50}'.format(' Finding any rois that span an entire contig '))
    start_time = time.time()
    end_df = roi_df[((roi_df['start_pos'].between(0,2000)) | (roi_df['end_pos'].between(roi_df['contig_len'] - 2000, roi_df['contig_len'] + 1)))].copy() # pull out rois into a seperate df that would have failed later tests unfairly due to being close to end of contig
    roi_df = roi_df[~(((roi_df['start_pos'].between(0,2000)) | (roi_df['end_pos'].between(roi_df['contig_len'] - 2000, roi_df['contig_len'] + 1))))].copy()
    for index, row in end_df.iterrows():
        if (0 <= row['start_pos'] <= 2000) and (row['contig_len'] - 2000 <= row['end_pos'] <= row['contig_len']):
            end_df.loc[index, 'end_info'] = 'whole'
        elif (0 <= row['start_pos'] <= 2000) and (not row['contig_len'] - 2000 <= row['end_pos'] <= row['contig_len']):
            end_df.loc[index, 'end_info'] = 'left'
        elif (not 0 <= row['start_pos'] <= 2000) and (row['contig_len'] - 2000 <= row['end_pos'] <= row['contig_len']):
            end_df.loc[index, 'end_info'] = 'right'
    roi_df['contig_split'] = np.nan
    roi_df['end_info'] = np.nan
    if len(roi_df) < 1:
        roi_df = 'exit'
    if len(end_df) < 1:
        end_df = 'exit'
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return roi_df, end_df

def process_clip_sam(output_folder,roi_df):
    print('\n{:#^50}'.format(' Processing clipped sam '))
    start_time = time.time()
    samfile = pysam.AlignmentFile(output_folder + '/temp_minimap_clipped_sorted.bam', 'rb')
    df_dict = {'qname':[], 'rname':[], 'pos':[], 'cigar':[],'tlen':[],'seq':[],'flag':[],'pnext':[]}
    for index, row in roi_df.iterrows():
        contig = row['contig']
        start = row['start_pos'] - 200
        if start < 0:
            start = 0
        end = row['end_pos'] + 200
        if end > row['contig_len']:
            end = row['contig_len']
        for reads in samfile.fetch(contig, start, end):
            df_dict['qname'].append(reads.query_name)
            df_dict['rname'].append(contig)
            df_dict['pos'].append(reads.reference_start)
            df_dict['cigar'].append(reads.cigarstring)
            df_dict['tlen'].append(reads.template_length)
            df_dict['seq'].append(reads.query_sequence)
            df_dict['flag'].append(reads.flag)
            df_dict['pnext'].append(reads.next_reference_start)
    df = pd.DataFrame.from_dict(df_dict)
    df = df.drop_duplicates(subset = ['qname', 'rname', 'pos', 'cigar', 'tlen', 'seq', 'flag', 'pnext'], keep='first')
    df['tlen'] = df['tlen'].astype(int)
    df['pos'] = df['pos'].astype(int)
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return df

def check_read_lengths(row, index, clip_df, side_name, names, side):
    df_dict = {'qname':[], 'rname':[], 'pos':[], 'cigar':[],'tlen':[],'seq':[],'flag':[],'pnext':[], 'roi':[], 'seq_len':[]}
    start_list = np.core.defchararray.add('{}~{}__'.format(row['roi'], str(row['start_pos'])), np.arange(1, row['start_count']+1, 1).astype(int).astype(str))
    if side == 'start':
        for read in start_list:
            try:
                start_reads = names.find(read)
                for x in start_reads:
                    df_dict['qname'].append(x.query_name)
                    df_dict['rname'].append(x.reference_name)
                    df_dict['pos'].append(x.reference_start)
                    df_dict['cigar'].append(x.cigarstring)
                    df_dict['tlen'].append(x.template_length)
                    df_dict['seq'].append(x.query_sequence)
                    df_dict['flag'].append(x.flag)
                    df_dict['pnext'].append(x.next_reference_start)
                    df_dict['roi'].append(row['roi'])
                    df_dict['seq_len'].append(len(x.query_sequence))
            except:
                pass
        start_df = pd.DataFrame.from_dict(df_dict)
        start_df = start_df.drop_duplicates()
        if len(start_df[start_df['seq_len'] >= 35]) < 10: ##### this is to make sure that any reads that didnt have clipped sections that were no long enough to be mapped are not unfairly excluded from the search
            return 10
    elif side == 'end':
        end_list = np.core.defchararray.add('{}~{}__'.format(row['roi'], str(row['end_pos'])), np.arange(row['start_count'], row['start_count'] + row['end_count']+1, 1).astype(int).astype(str))
        for read in end_list:
            try:
                end_reads = names.find(read)
                for x in end_reads:
                    df_dict['qname'].append(x.query_name)
                    df_dict['rname'].append(x.reference_name)
                    df_dict['pos'].append(x.reference_start)
                    df_dict['cigar'].append(x.cigarstring)
                    df_dict['tlen'].append(x.template_length)
                    df_dict['seq'].append(x.query_sequence)
                    df_dict['flag'].append(x.flag)
                    df_dict['pnext'].append(x.next_reference_start)
                    df_dict['roi'].append(row['roi'])
                    df_dict['seq_len'].append(len(x.query_sequence))
            except:
                pass
        end_df = pd.DataFrame.from_dict(df_dict)
        end_df = end_df.drop_duplicates()
        if len(end_df[end_df['seq_len'] >= 35]) < 10: ##### this is to make sure that any reads that didnt have clipped sections that were no long enough to be mapped are not unfairly excluded from the search
            return 10


def get_clip_pos(clip_df,roi_df, output_folder):
    print('\n{:#^50}'.format(' Getting pos for clipped ends (contained rois) '))
    start_time = time.time()
    samfile = pysam.AlignmentFile(output_folder + '/temp_minimap_clipped_sorted.bam', 'rb')
    names = pysam.IndexedReads(samfile)
    names.build()
    roi_df['start-end_clip_count'] = np.nan
    roi_df['end-start_clip_count'] = np.nan
    roi_df['total_clip_count'] = np.nan
    clip_df2 = clip_df[(clip_df['cigar'].str.count('M') == 1) & (clip_df['cigar'].str.count('S') == 0) & (clip_df['cigar'].str.count('I') == 0) & (clip_df['cigar'].str.count('X') == 0) & (clip_df['cigar'].str.count('D') == 0)]
    for index, row in roi_df.iterrows():
        start_name = row['roi'] + '~' + str(row['start_pos'])
        end_name = row['roi'] + '~' + str(row['end_pos'])
        start_mapped_df = clip_df2[(clip_df2['qname'].str.split('__').str[0] == start_name) & (clip_df2['pos'].between(row['end_pos'] - 200, row['end_pos'] + 200))]
        end_mapped_df = clip_df2[(clip_df2['qname'].str.split('__').str[0] == end_name) & (clip_df2['pos'].between(row['start_pos'] - 200, row['start_pos'] + 200))]
        if (len(start_mapped_df) >= 10) and (len(end_mapped_df) >= 10):
            roi_df.loc[index, 'start-end_clip_count'] = len(end_mapped_df)
            roi_df.loc[index, 'end-start_clip_count'] = len(start_mapped_df)
            roi_df.loc[index, 'total_clip_count'] = len(end_mapped_df) + len(start_mapped_df) + row['start_count'] + row['end_count'] # added these last 2 to test something
        else:
            if (len(start_mapped_df) < 10) and (len(end_mapped_df) >= 10):
                roi_df.loc[index, 'start-end_clip_count'] = len(end_mapped_df)
                roi_df.loc[index, 'end-start_clip_count'] = check_read_lengths(row, index, clip_df, start_name, names, 'start')
                roi_df.loc[index, 'total_clip_count'] = 10 + len(start_mapped_df) + row['start_count'] + row['end_count'] # added these last 2 to test something
            elif (len(end_mapped_df) < 10) and (len(start_mapped_df) >= 10):
                roi_df.loc[index, 'end-start_clip_count'] = len(start_mapped_df)
                roi_df.loc[index, 'start-end_clip_count'] = check_read_lengths(row, index, clip_df, end_name, names, 'end')
                roi_df.loc[index, 'total_clip_count'] = len(end_mapped_df) + 10 + row['start_count'] + row['end_count'] # added these last 2 to test something
    roi_df = roi_df[(roi_df['start-end_clip_count'].notna()) & (roi_df['end-start_clip_count'].notna())]
    if len(roi_df) < 1:
        return 'exit'
    roi_df = roi_df.sort_values(by=['roi','total_clip_count'], ascending=False)
    df_list = []
    for i in roi_df['roi'].unique():
        df_list.append(roi_df[roi_df['roi'] == i].iloc[[0]])
    roi_df = pd.concat(df_list)
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return roi_df

def process_clip_sam_end_rois(output_folder,roi_df):
    print('\n{:#^50}'.format(' Processing clipped sam for ends '))
    start_time = time.time()
    samfile = pysam.AlignmentFile(output_folder + '/temp_minimap_clipped_sorted.bam', 'rb')
    names = pysam.IndexedReads(samfile)
    names.build()
    df_dict = {'qname':[], 'rname':[], 'pos':[], 'cigar':[],'tlen':[],'seq':[],'flag':[],'pnext':[], 'roi':[]}
    for index, row in roi_df.iterrows():
        start_list = np.core.defchararray.add('{}~{}__'.format(row['roi'], str(row['start_pos'])), np.arange(1, row['start_count']+1, 1).astype(int).astype(str))
        for read in start_list:
            try:
                start_reads = names.find(read)
                for x in start_reads:
                    df_dict['qname'].append(x.query_name)
                    df_dict['rname'].append(x.reference_name)
                    df_dict['pos'].append(x.reference_start)
                    df_dict['cigar'].append(x.cigarstring)
                    df_dict['tlen'].append(x.template_length)
                    df_dict['seq'].append(x.query_sequence)
                    df_dict['flag'].append(x.flag)
                    df_dict['pnext'].append(x.next_reference_start)
                    df_dict['roi'].append(row['roi'])
            except:
                pass
        end_list = np.core.defchararray.add('{}~{}__'.format(row['roi'], str(row['end_pos'])), np.arange(row['start_count'], row['start_count'] + row['end_count']+1, 1).astype(int).astype(str))
        for read in end_list:
            try:
                end_reads = names.find(read)
                for x in end_reads:
                    df_dict['qname'].append(x.query_name)
                    df_dict['rname'].append(x.reference_name)
                    df_dict['pos'].append(x.reference_start)
                    df_dict['cigar'].append(x.cigarstring)
                    df_dict['tlen'].append(x.template_length)
                    df_dict['seq'].append(x.query_sequence)
                    df_dict['flag'].append(x.flag)
                    df_dict['pnext'].append(x.next_reference_start)
                    df_dict['roi'].append(row['roi'])
            except:
                pass
    df = pd.DataFrame.from_dict(df_dict)
    df = df.drop_duplicates(subset = ['qname', 'rname', 'pos', 'cigar', 'tlen', 'seq', 'flag', 'pnext'], keep='first')
    df = df.dropna(subset=['rname'])
    df['tlen'] = df['tlen'].astype(int)
    df['pos'] = df['pos'].astype(int)
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return df

def get_clip_pos_end_rois(clip_df,roi_df,seq):
    print('\n{:#^50}'.format('  Getting pos for clipped ends (end rois) '))
    start_time = time.time()
    clip_df = clip_df[(clip_df['cigar'].str.count('M') == 1) & (clip_df['cigar'].str.count('S') == 0) & (clip_df['cigar'].str.count('I') == 0) & (clip_df['cigar'].str.count('X') == 0) & (clip_df['cigar'].str.count('D') == 0)]
    sub_df = clip_df['qname'].str.split('__').str[0].value_counts().rename_axis('coords').reset_index(name='counts')
    df_list = []
    for index, row in clip_df.iterrows():
        clip_df.loc[index,'slen'] = len(row['seq'])
    for index, row in roi_df.iterrows():
        sub_list = []
        start_name = row['roi'] + '~' + str(row['start_pos'])
        end_name = row['roi'] + '~' + str(row['end_pos'])
        start_df = clip_df[clip_df['qname'].str.split('__').str[0] == start_name].copy()
        end_df = clip_df[clip_df['qname'].str.split('__').str[0] == end_name].copy()
        start_df['query'] = start_df['qname'].str.split('__').str[0].copy()
        start_df = start_df.groupby(['query', 'rname'])['qname'].count()
        start_df = start_df[start_df >= 10]
        end_df['query'] = end_df['qname'].str.split('__').str[0].copy()
        end_df = end_df.groupby(['query', 'rname'])['qname'].count()
        end_df = end_df[end_df >= 10]
        if (len(start_df) > 0) and (len(end_df) > 0):
            start_df = start_df.reset_index(name='counts')
            for i in start_df['rname'].unique():
                sub_df = clip_df[(clip_df['qname'].str.split('__').str[0] == start_name) & (clip_df['rname'] == i)].copy()
                start_contig_length = len(seq[i].seq)
                sub_df2 = sub_df.groupby(pd.cut(sub_df['pos'], np.arange(-1,start_contig_length + 201 ,200)))['qname'].count()
                sub_df2 = sub_df2[sub_df2 >= 10]
                sub_df2 = sub_df2.reset_index(name='counts')
                sub_df2['start_query_pos'] = start_df['query'].str.split('~').str[-1].astype(int)
                sub_df2['query'] = start_df['query'].str.split('~').str[:4].str.join('~')
                sub_df2 = sub_df2.sort_values(by=['counts'], ascending=False).iloc[[0]]
                for ind, r in sub_df2.iterrows():
                    if not isinstance(r['pos'],float):
                        sub_df2.loc[ind, 'left'] = r['pos'].left
                        sub_df2.loc[ind, 'right'] = r['pos'].right
                    else:
                        sub_df2.loc[ind, 'left'] = np.nan
                        sub_df2.loc[ind, 'right'] = np.nan
                sub_df = sub_df[sub_df['pos'].between(sub_df2['left'].iloc[0],sub_df2['right'].iloc[0])]
                roi_df.loc[index, 'right_end_pos'] = np.median(sub_df['pos'])
                roi_df.loc[index, 'right_end_count'] = sub_df2['counts'].iloc[0]
                roi_df.loc[index, 'right_end_rname'] = i
            end_df = end_df.reset_index(name='counts')
            for i in end_df['rname'].unique():
                sub_df = clip_df[(clip_df['qname'].str.split('__').str[0] == end_name) & (clip_df['rname'] == i)].copy()
                end_contig_length = len(seq[i].seq)
                sub_df2 = sub_df.groupby(pd.cut(sub_df['pos'], np.arange(-1,end_contig_length + 201 ,200)))['qname'].count()
                sub_df2 = sub_df2[sub_df2 >= 10]
                sub_df2 = sub_df2.reset_index(name='counts')
                sub_df2['end_query_pos'] = end_df['query'].str.split('~').str[-1].astype(int)
                sub_df2['query'] = end_df['query'].str.split('~').str[:4].str.join('~')
                sub_df2 = sub_df2.sort_values(by=['counts'], ascending=False).iloc[[0]]
                for ind, r in sub_df2.iterrows():
                    sub_df2.loc[ind, 'left'] = r['pos'].left
                    sub_df2.loc[ind, 'right'] = r['pos'].right
                sub_df = sub_df[sub_df['pos'].between(sub_df2['left'].iloc[0],sub_df2['right'].iloc[0])]
                roi_df.loc[index, 'left_end_pos'] = np.median(sub_df['pos'])
                roi_df.loc[index, 'left_end_count'] = sub_df2['counts'].iloc[0]
                roi_df.loc[index, 'left_end_rname'] = i
        else:
            roi_df.loc[index, 'right_end_pos'] = np.nan
            roi_df.loc[index, 'right_end_count'] = np.nan
            roi_df.loc[index, 'right_end_rname'] = np.nan
            roi_df.loc[index, 'left_end_pos'] = np.nan
            roi_df.loc[index, 'left_end_count'] = np.nan
            roi_df.loc[index, 'left_end_rname'] = np.nan
    roi_df = roi_df.dropna()
    if len(roi_df) < 1:
        return 'exit'
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return roi_df

def filter_by_end_status(roi_df):
    print('\n{:#^50}'.format(' Filtering end rois '))
    start_time = time.time()
    for index, row in roi_df.iterrows():
        if (row['end_info'] == 'whole') and ((row['contig'] == row['right_end_rname']) or (row['contig'] == row['left_end_rname'])):
            roi_df.drop(index, inplace=True)
        elif (row['end_info'] == 'left') and (row['contig'] == row['left_end_rname']):
            roi_df.drop(index, inplace=True)
        elif (row['end_info'] == 'right') and (row['contig'] == row['right_end_rname']):
            roi_df.drop(index, inplace=True)
    if len(roi_df) < 1:
        return 'exit'
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return roi_df

def fix_table_whole(index, row, roi_df, seq):
    left_contig_length = len(seq[row['left_end_rname']].seq)
    closest_left = min([0,left_contig_length], key=lambda x:abs(x-row['left_end_pos']))
    if closest_left == 0:
        left_start = row['left_end_pos']
        left_end = 0
    else:
        left_start = left_contig_length
        left_end = row['left_end_pos']
    right_contig_length = len(seq[row['right_end_rname']].seq)
    closest_right = min([0,right_contig_length], key=lambda x:abs(x-row['right_end_pos']))
    if closest_right == 0:
        right_start = 0
        right_end = row['right_end_pos']
    else:
        right_start = row['right_end_pos']
        right_end = righ_contig_length
    roi_df.loc[index, 'contig_split'] = '{}({}..{}) -> {}({}..{}) -> {}({}..{})'.format(row['left_end_rname'],left_start,left_end,
                                                                                     row['contig'],row['start_pos'],row['end_pos'],
                                                                                     row['right_end_rname'],right_start,right_end)
    return roi_df

def fix_table_left(index, row, roi_df, seq):
    left_contig_length = len(seq[row['left_end_rname']].seq)
    closest_left = min([0,left_contig_length], key=lambda x:abs(x-row['left_end_pos']))
    if closest_left == 0:
        left_start = row['left_end_pos']
        left_end = 0
    else:
        left_start = left_contig_length
        left_end = row['left_end_pos']
    roi_df.loc[index, 'contig_split'] = '{}({}..{}) -> {}({}..{})'.format(row['left_end_rname'],left_start,left_end,
                                                                                     row['contig'],row['start_pos'],row['end_pos'])
    return roi_df

def fix_table_right(index, row, roi_df, seq):
    right_contig_length = len(seq[row['right_end_rname']].seq)
    closest_right = min([0,right_contig_length], key=lambda x:abs(x-row['right_end_pos']))
    if closest_right == 0:
        right_start = row['right_end_pos']
        right_end = 0
    else:
        right_start = right_contig_length
        right_end = row['right_end_pos']
    rroi_df.loc[index, 'contig_split'] = '{}({}..{}) -> {}({}..{})'.format(row['contig'],row['start_pos'],row['end_pos'],
                                                                                     row['right_end_rname'],right_start,right_end)
    return roi_df


def reformat_end_roi_tables(roi_df, seq):
    print('\n{:#^50}'.format(' Reformatting end roi table '))
    start_time = time.time()
    roi_df['start_pos'], roi_df['end_pos'], roi_df['left_end_pos'], roi_df['right_end_pos'] = roi_df['start_pos'].astype(int), roi_df['end_pos'].astype(int), roi_df['left_end_pos'].astype(int), roi_df['right_end_pos'].astype(int)
    for index, row in roi_df.iterrows():
        if row['end_info'] == 'whole':
            roi_df = fix_table_whole(index, row, roi_df,seq)
        elif row['end_info'] == 'left':
            roi_df = fix_table_left(index, row, roi_df,seq)
        elif row['end_info'] == 'right':
            roi_df = fix_table_right(index, row, roi_df,seq)
    roi_df = roi_df[['roi', 'start_pos', 'start_count', 'end_pos', 'end_count', 'length', 'contig', 'contig_len', 'end_info', 'contig_split', 'med_z']]
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return roi_df

def get_longest_len_z(z, cutoff):
    x = {'start': [], 'end': []}
    above_threshold = np.diff(np.array(z) < cutoff, prepend=False) # less than as this then will flag start and end regions as where the go below and come above
    rise_above_thres = np.argwhere(above_threshold)[::2,0]
    sink_below_thres = np.argwhere(above_threshold)[1::2,0]
    x['start'] = rise_above_thres
    x['end'] = sink_below_thres
    lengths = []
    if len(x['start']) > 0:
        for j in np.arange(0,len(x['start']) - 1 ):
            lengths.append(x['end'][j] - x['start'][j])
        if len(lengths) > 0:
            longest_below_z = max(lengths)
        else:
            longest_below_z = 0
    else:
        longest_below_z = 0
    return longest_below_z


def multiprocessing_calc_roi_z_medians(roi_df,depths,median,mad,median_z_cutoff):
    for i in roi_df['contig'].unique():
        z = [((0.6745*(x - median))/mad) for x in depths[i]]
        for index, row in roi_df[(roi_df['contig'] == i)].iterrows():
            med_z = z[row['start_pos']:row['end_pos']]
            med_z = np.median(med_z)
            med_z_l = z[(row['start_pos'] - 10000):row['start_pos']]
            longest_below_z = get_longest_len_z(z[row['start_pos']:row['end_pos']], (med_z/2))
            if len(med_z_l) > 0:
                med_z_l = np.median(med_z_l)
            else:
                med_z_l = np.nan
            med_z_r = z[row['end_pos']:(row['end_pos'] + 10000)]
            if len(med_z_r) > 0:
                med_z_r = np.median(med_z_r)
            else:
                med_z_r = np.nan
            roi_df.loc[index, 'med_z'] = med_z
            roi_df.loc[index, 'med_z_l'] = med_z_l
            roi_df.loc[index, 'med_z_r'] = med_z_r
            roi_df.loc[index,'longest_below_z'] = longest_below_z
    return roi_df

def calc_roi_Z_medians(roi_df,depths,median,mad,median_z_cutoff,threads):
    print('\n{:#^50}'.format(' Calculating median Z for each roi '))
    start_time = time.time()
    roi_df = roi_df.reset_index(drop=True)
    roi_df['start_pos'] = roi_df['start_pos'].astype('int')
    roi_df['end_pos'] = roi_df['end_pos'].astype('int')
    chunks, pool = hZ.misc_functions.split_df_for_multiproc(roi_df,threads)
    roi_df = pd.concat(pool.map(partial(multiprocessing_calc_roi_z_medians, depths=depths, median=median, mad=mad, median_z_cutoff=median_z_cutoff), chunks))
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return roi_df

def filter_Z_medians(roi_df,median_z_cutoff):
    print('\n{:#^50}'.format(' Filtering rois by Z-score '))
    start_time = time.time()
    both_ends_df = roi_df[(roi_df['end_info'] == 'whole') & (roi_df['longest_below_z'] < 7500)]
    left_df = roi_df[(roi_df['end_info'] == 'left') & (roi_df['longest_below_z'] < 7500) &
                        ((roi_df['med_z'] >= roi_df['med_z']/2) & (roi_df['med_z_r'] < roi_df['med_z']/2) & (roi_df['med_z_l'].isna()) |
                        (roi_df['med_z'] >= roi_df['med_z']/2) & (roi_df['med_z_r'] < roi_df['med_z']/2) & (roi_df['med_z_l'] < roi_df['med_z']/2))]
    right_df = roi_df[(roi_df['end_info'] == 'right') & (roi_df['longest_below_z'] < 7500) &
                        ((roi_df['med_z'] >= roi_df['med_z']/2) & (roi_df['med_z_l'] < roi_df['med_z']/2) & (roi_df['med_z_r'].isna()) |
                        (roi_df['med_z'] >= roi_df['med_z']/2) & (roi_df['med_z_r'] < roi_df['med_z']/2) & (roi_df['med_z_l'] < roi_df['med_z']/2))]
    contained_df = roi_df[(roi_df['med_z'] >= roi_df['med_z']/2) & (roi_df['longest_below_z'] < 7500) & (roi_df['med_z_l'] < roi_df['med_z']/2) & (roi_df['med_z_r'] < roi_df['med_z']/2) & (roi_df['end_info'].isna())]
    roi_df = pd.concat([both_ends_df,left_df, right_df, contained_df])
    roi_df = roi_df.reset_index(drop=True)
    if len(roi_df) < 1:
        return 'exit'
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return roi_df

def check_for_circular_roi(sam_df,roi_df):
    print('\n{:#^50}'.format(' Checking for plasmid phages '))
    start_time = time.time()
    overlaps_dict = {}
    roi_df['circular'] = False
    for index,row in roi_df[~pd.isna(roi_df['end_info'])].iterrows():
        start = row['start_pos']
        end = row['end_pos']
        contig_len = row['contig_len']
        contig = row['roi'].split('~')[0]
        read_list = []
        if start < 1000:
            read_list = list(sam_df['qname'][(sam_df['rname'] == contig) & (sam_df['pos'].between(0,1000))])
        elif end > contig_len - 1000:
            read_list = list(sam_df['qname'][(sam_df['rname'] == contig) & (sam_df['pos'].between(contig_len - 1000,contig_len))])
        sub_df = sam_df[(sam_df['qname'].isin(read_list)) & (sam_df['rname'] == contig)]
        sub_df = sub_df[((sub_df['pos'].between(0,1000)) & (sub_df['pnext'].between(contig_len - 1000,contig_len))) | ((sub_df['pnext'].between(0,1000)) & (sub_df['pos'].between(contig_len - 1000,contig_len)))]
        sub_df = sub_df[sub_df.duplicated(subset = ['qname','rname'], keep=False)]
        min = sub_df['pos'][sub_df['pos'] <= 1000].max()
        max = sub_df['pos'][sub_df['pos'] >= contig_len-1000].min()
        overlaps = roi_df['roi'][(roi_df['roi'].str.split('~').str[0] == contig) & ((roi_df['start_pos'].between(0,min)) | (roi_df['end_pos'].between(max,contig_len)))].unique()
        if len(overlaps) >= 1:
            overlaps_dict[row['roi']] = overlaps
    overlaps = []
    for i in overlaps_dict:
        overlaps.append([k for (k, v) in overlaps_dict.items() if i in v])
    unique_overlaps = [list(x) for x in set(tuple(x) for x in overlaps)]
    counter = 1000
    for overlaps in unique_overlaps:
        starts = []
        ends = []
        for i in overlaps:
            starts.append(roi_df['start_pos'][roi_df['roi'] == i].iloc[0])
            ends.append(roi_df['end_pos'][roi_df['roi'] == i].iloc[0])
            contig_len = roi_df['contig_len'][roi_df['roi'] == i].iloc[0]
            roi_df = roi_df[roi_df['roi'] != i]
        contig = i.split('~')[0]
        start = np.max(starts)
        end = np.min(ends)
        roi_1 = contig + '~roi_' + str(counter) + '_c1'
        roi_2 = contig + '~roi_' + str(counter) + '_c2'
        df2 = {'roi': [roi_1,roi_2], 'start_pos': [start,end], 'start_count': [np.nan,np.nan],'end_pos': [end,start],'end_count': [np.nan,np.nan], 'contig_len': [contig_len,contig_len], 'contig': [contig,contig], 'circular':[True,True]}
        df2 = pd.DataFrame.from_dict(df2)
        roi_df = pd.concat([roi_df,df2])
        counter = counter - 1
    circular_df = roi_df[roi_df['circular'] == True]
    roi_df = roi_df[roi_df['circular'] == False]
    if len(roi_df) < 1:
        roi_df = 'exit'
    if len(circular_df) < 1:
        circular_df = 'exit'
    else:
        for index, row in circular_df.iterrows():
            circular_df.loc[index, 'contig_split'] = np.nan
    end_time = '{:.3f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return circular_df, roi_df

def extract_roi_orfs(orf_df, orf_aa, roi_df, orf_dna, output_folder,min_orfs):
    print('\n{:#^50}'.format(' Getting orf sequences for rois '))
    start_time = time.time()
    roi_aa = {}
    roi_dna = {}
    all_roi_orfs = []
    orf_no = []
    for index, row in roi_df.iterrows():
        contig = row['roi'].split('~')[0]
        counter = 1
        if (row['circular'] == False) | ((row['circular'] == True) & (row['roi'].split('_')[-1] == 'c1')):
            orfs = orf_df[(orf_df['contig'] == contig) & (orf_df['start'] >= row['start_pos']) & (orf_df['end'] <= row['end_pos'])].copy()
        elif (row['circular'] == True) & (row['roi'].split('_')[-1] == 'c2'):
            orfs_1 = orf_df[(orf_df['contig'] == contig) & (orf_df['start'] >= row['start_pos']) & (orf_df['end'] <= row['contig_len'])].copy()
            orfs_2 = orf_df[(orf_df['contig'] == contig) & (orf_df['start'] >= 0) & (orf_df['end'] <= row['end_pos'])].copy()
            orfs = pd.concat([orfs_1,orfs_2])

        if len(orfs) < min_orfs:
            roi_df = roi_df[roi_df['roi'] != row['roi']]
        else:
            for i,r in orfs.iterrows():
                orfs.loc[i, 'new_orf'] = 'orf_' + str(counter)
                counter = counter + 1
            orf_no.append(len(orfs))
            roi_aa[row['roi']] = []
            roi_dna[row['roi']] = []
            for i in list(orfs['orf']):
                aa = orf_aa[contig + '_' + i]
                aa.id = row['roi'] + '~' + orfs['new_orf'][orfs['orf'] == i].iloc[0]
                roi_aa[row['roi']].append(aa)
                dna = orf_dna[contig + '_' + i]
                dna.id = row['roi'] + '~' + orfs['new_orf'][orfs['orf'] == i].iloc[0]
                roi_dna[row['roi']].append(dna)
                all_roi_orfs.append(aa)
    if len(roi_df) < 1:
        return 'exit', 'exit', 'exit'
    else:
        roi_df['orf_count'] = orf_no.copy()
        roi_df = roi_df.reset_index(drop = True)
        SeqIO.write(all_roi_orfs, output_folder + '/temp_aa.fasta', 'fasta')
        end_time = '{:.3f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return roi_df, roi_aa, roi_dna

def hmm_scan(output_folder,threads,db_path):
    print('\n{:#^50}'.format(' Screening orfs vs pVOGs HMM files '))
    start_time = time.time()
    aa_file = output_folder + '/temp_aa.fasta'
    hmm_out = output_folder + '/temp_hmmout.txt'
    hmm_db_path = db_path + '/combined_AllvogHMMprofiles.hmm'
    hmm_table_path = db_path + '/combined_Allvogtables.txt'
    subprocess.run(['hmmscan',
                    '--noali',
                    '--notextw',
                    '--cpu',
                    str(threads),
                    '-E',
                    '0.00001',
                    '--tblout',
                    hmm_out,
                    hmm_db_path,
                    aa_file], stdout=subprocess.DEVNULL)
    df = pd.read_csv(hmm_out, comment = '#', delim_whitespace=True, usecols = [0,2,4], names = ('vog','orf','e'))
    df = df.iloc[df.groupby(['orf'])['e'].idxmin()].reset_index(drop=True)
    df['orf_no'] = df['orf'].str.split('~').str[-1]
    df = df.sort_values(by='orf_no')
    df2 = pd.read_csv(hmm_table_path, index_col=0)
    df = pd.merge(df,df2[['vog','type','description']], on='vog', how = 'left')
    df.to_csv(output_folder + '/temp_hmms.tab', sep = '\t', index = False)
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return df

def calc_pVOGs_frac(hmm_df,roi_df,frac_pvog):
    print('\n{:#^50}'.format(' Calculating fraction of orfs hit by pVOGs '))
    start_time = time.time()
    hmm_df['prophage'] = hmm_df['orf'].str.split('~').str[:-1].str.join('~')
    hmm_counts = hmm_df['prophage'].value_counts().reset_index()
    hmm_counts.columns = ['roi','hmm']
    frac_df = pd.merge(roi_df,hmm_counts,on='roi', how = 'left')
    frac_df['frac_pvog'] = frac_df['hmm']/frac_df['orf_count']
    df = pd.merge(roi_df,frac_df[['roi','frac_pvog']], on ='roi', how = 'left')
    df['frac_pvog'] = df['frac_pvog'].copy().fillna(0)
    df = df[df['frac_pvog'] >= frac_pvog]
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    if len(df) < 1:
        return 'exit'
    return df

def get_att(roi_df,seq,output_folder):
    print('\n{:#^50}'.format(' Get att sites '))
    start_time = time.time()
    atts = {'roi': [], 'attL_pos': [], 'attR_pos': [], 'attL_seq': [], 'attR_seq': []}
    left_list = []
    right_list = []
    roi_df['start_pos'] = roi_df['start_pos'].astype(int)
    roi_df['end_pos'] = roi_df['end_pos'].astype(int)
    for index, row in roi_df.iterrows():
        if (row['circular'] == True):
            roi_df.loc[index, 'attL_seq'] = np.nan
            roi_df.loc[index, 'attR_seq'] = np.nan
        if not pd.isna(row['contig_split']):
            roi_df.loc[index, 'attL_seq'] = np.nan
            roi_df.loc[index, 'attR_seq'] = np.nan
        else:
            left = seq[row['contig']][row['start_pos'] - 100 :row['start_pos'] + 100]
            right = seq[row['contig']][row['end_pos'] - 100 :row['end_pos'] + 100]
            left.id = row['roi']
            right.id = row['roi']
            left_list.append(left)
            right_list.append(right)
    if len(left_list) > 0 and len(right_list) > 0:
        SeqIO.write(left_list, output_folder + '/temp_left.fasta', 'fasta')
        SeqIO.write(right_list, output_folder + '/temp_right.fasta', 'fasta')
        p = subprocess.run(['blastn',
                                '-query',
                                output_folder + '/temp_left.fasta',
                                '-subject',
                                output_folder + '/temp_right.fasta',
                                '-evalue',
                                '10000',
                                '-task',
                                'blastn-short',
                                '-outfmt',
                                '6 qseqid qstart qend sseqid sstart send evalue qseq sseq length'], stdout = subprocess.PIPE)
        s=str(p.stdout,'utf-8')
        blast = pd.read_csv(StringIO(s),sep="\t",names=('qname','qstart','qend','sname','sstart','send','e','qseq','sseq','length'))
        blast = blast.astype({'e': float})
        blast = blast[blast['length'] > 11]
    else:
        column_names = ['qname','qstart','qend','sname','sstart','send','e','qseq','sseq','length']
        blast = pd.DataFrame(columns = column_names)
    for index, row in roi_df.iterrows():
        sub_df = blast[(blast['qname'] == row['roi']) & (blast['sname'] == row['roi'])]
        if len(sub_df) > 0:
            sub_df = sub_df.sort_values(by = ['e'])
            roi_df.loc[index, 'attL_seq'] = sub_df['qseq'].iloc[0]
            roi_df.loc[index, 'attR_seq'] = sub_df['sseq'].iloc[0]
        else:
            roi_df.loc[index, 'attL_seq'] = np.nan
            roi_df.loc[index, 'attR_seq'] = np.nan
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return roi_df

def quick_filter(roi_df):
    print('\n{:#^50}'.format(' Doing quick filter as no_extra == True '))
    start_time = time.time()
    df_list = []
    for i in roi_df['roi'].unique():
        sub_df = roi_df[roi_df['roi'] == i].copy()
        for index, row in sub_df.iterrows():
            if row['circular'] != True:
                sub_df.loc[index,'total_clips'] = int(row['start_count']) + int(row['end_count'])
            else:
                sub_df.loc[index,'total_clips'] = np.nan
        sub_df = sub_df.sort_values(by=['total_clips'],ascending=False).iloc[[0]]
        df_list.append(sub_df)
    roi_df = pd.concat(df_list)
    roi_df = roi_df.reset_index(drop=True)
    roi_df['contig_split'] = np.nan
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return roi_df

def keep_only_x_best(roi_df,keep_threshold):
    print('\n{:#^50}'.format(' Keeping only "best" roi possibilities '))
    start_time = time.time()
    df_list = []
    for i in roi_df['roi'].unique():
        sub_df = roi_df[roi_df['roi'] == i]
        sub_df = sub_df.reset_index(drop=True)
        for index, row in sub_df.iterrows():
            sub_df.loc[index,'total_clips'] = int(row['start_count']) + int(row['end_count'])
        sub_df = sub_df.sort_values(by=['total_clips'],ascending=False).iloc[0:keep_threshold,:]
        df_list.append(sub_df)
    roi_df = pd.concat(df_list)
    roi_df = roi_df.reset_index(drop=True)
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return roi_df

def split_df_for_multiproc_phrogs(lst,threads):
    chunks = round(len(lst)/threads)
    chunks = [lst[i:i + chunks] for i in range(0, len(lst), chunks)]
    pool = multiprocessing.Pool(processes=threads)
    return chunks, pool

def hhblits_func(chunk, path, out):
    for i in chunk:
        name = i.id
        i.id = '~'.join(i.id.split('~')[-2:])
        fasta_seq = ">{}\n{}".format(i.id, i.seq)
        subprocess.run(['hhblits', 
                        '-i', 'stdin',
                        '-blasttab', out + '/temp_' + i.id + '.tab',
                        '-n', '1',
                        '-cpu', '1',
                        '-d', path + '/phrogs', 
                        '-v', '0'], input = fasta_seq.encode('utf-8'), stdout=subprocess.DEVNULL)
        df = pd.read_csv(out + '/temp_' + i.id + '.tab', sep='\t',  usecols = [0,1,10], names = ['orf','phrog','e'])
        df = df.iloc[df.groupby(['orf'])['e'].idxmin()].reset_index(drop=True)
        df.loc[0, 'orf'] = name 
        df['orf_no'] = df['orf'].str.split('~').str[-1]
        df.to_csv(out + '/temp_' + i.id + '.tab', sep='\t', index=False, header=False)

def hhblits_phrogs(output_folder, threads, db_path):
    print('\n{:#^50}'.format(' Screening orfs vs pVOGs HMM files '))
    start_time = time.time()
    aa_file = output_folder + '/temp_aa.fasta'
    hhblits_out = output_folder + '/temp_hmmout.txt'
    hhblits_db_path = db_path + 'phrogs_hhsuite_db'
    hhblits_table_path = db_path + '/phrogs_table_almostfinal_plusGO_wNA_utf8.tsv'
    aa_list = []
    for i in SeqIO.parse(aa_file, 'fasta'):
        aa_list.append(i)
    chunks, pool = split_df_for_multiproc_phrogs(aa_list,threads)
    pool.map(partial(hhblits_func, path = hhblits_db_path, out = output_folder), chunks)
    os.system("cat {}/temp_*orf*.tab > {}/temp_hmms.tab".format(output_folder,output_folder))
    df = pd.read_csv(output_folder + '/temp_hmms.tab', sep='\t',  names = ['orf','phrog','e','orf_no'])
    df = df.sort_values(by='orf_no')
    for index, row in df.iterrows():
        if row['e'] > 0.00001:
            df.loc[index, 'e'] = np.nan
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return df

def calc_phrogs_frac(hmm_df,roi_df,frac_pvog):
    print('\n{:#^50}'.format(' Calculating fraction of orfs hit by pVOGs '))
    start_time = time.time()
    hmm_df['prophage'] = hmm_df['orf'].str.split('~').str[:-1].str.join('~')
    hmm_hit = hmm_df[~np.isnan(hmm_df['e'])]
    hmm_counts = hmm_hit['prophage'].value_counts().reset_index()
    hmm_counts.columns = ['roi','hmm']
    frac_df = pd.merge(roi_df,hmm_counts,on='roi', how = 'left')
    frac_df['frac_pvog'] = frac_df['hmm']/frac_df['orf_count']
    df = pd.merge(roi_df,frac_df[['roi','frac_pvog']], on ='roi', how = 'left')
    df['frac_pvog'] = df['frac_pvog'].copy().fillna(0)
    df = df[df['frac_pvog'] >= frac_pvog]
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    if len(df) < 1:
        return 'exit'
    return df
