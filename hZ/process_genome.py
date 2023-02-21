import time
from Bio import SeqIO
import sys
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import pyrodigal


def process_fasta(fasta,output_folder):
    print('\n{:#^50}'.format(' Processing fasta '))
    start_time = time.time()
    raw_seq_dict = {}
    seq_dict = {}
    seq_list = []
    seq_temp_out = output_folder + '/temp_genome.fasta'
    for i in SeqIO.parse(fasta, 'fasta'):
        raw_seq_dict[i.id] = i
        if len(i.seq) > 10000:
            seq_dict[i.id] = i
            seq_list.append(i)
    if len(seq_dict) < 1:
        print('hafeZ: error: no contigs were of length that passed initial filtering, check fasta input is correct')
        return 'exit', 'exit', 'exit'
    elif len(seq_dict) == 1:
        multicontig = False
    elif len(seq_dict) > 1:
        multicontig = True
    SeqIO.write(seq_list, seq_temp_out, 'fasta')
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return seq_dict, multicontig, seq_temp_out, raw_seq_dict

def get_genome_length(raw_seq_dict):
    print('\n{:#^50}'.format(' Calculating genome length '))
    start_time = time.time()
    genome_length = 0
    for i in raw_seq_dict:
        genome_length = genome_length + len(raw_seq_dict[i].seq)
    end_time = '{:.5f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return genome_length

def get_roi_sequences(roi_df, seq, output_folder):
    print('\n{:#^50}'.format(' Extracting rois DNA sequences '))
    start_time = time.time()
    roi_list = []
    for index, row in roi_df.iterrows():
        if (row['circular'] == False) | ((row['circular'] == True) & (row['roi'].split('_')[-1] == 'c1')):
            start = int(row['start_pos'])
            end = int(row['end_pos'])
            roi = seq[row['roi'].split('~')[0]][start:end-1]
            roi.id = row['roi']
            roi.name = ""
            roi.description = ""
            roi_list.append(roi)
        elif (row['circular'] == True) & (row['roi'].split('_')[-1] == 'c2'):
            start_1 = int(row['start_pos'])
            end_1 = int(row['contig_len'])
            roi_1 = seq[row['roi'].split('~')[0]][start_1:end_1-1]
            start_2 = 0
            end_2 = int(row['end_pos'])
            roi_2 = seq[row['roi'].split('~')[0]][start_2:end_2-1]
            roi = SeqIO.SeqRecord(Seq("".join([str(seq_rec) for seq_rec in [roi_1.seq,roi_2.seq]])), id = row['roi'], name = '', description = '')
            roi_list.append(roi)
    SeqIO.write(roi_list, output_folder + "/temp_refined_roi.fasta", 'fasta')
    end_time = '{:.3f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return roi_list



def get_orfs(seq_dict, multicontig):
    print('\n{:#^50}'.format(' Getting ORFs '))
    start_time = time.time()
    orfs_aa = {}
    orfs_dna = {}
    orfs_df_list = []
    p = pyrodigal.OrfFinder(meta = multicontig)
    for i in seq_dict:
        if multicontig is False:
            p.train(str(seq_dict[i].seq))
        genes = p.find_genes(str(seq_dict[i].seq))
        orfs_db = {'contig': [], 'orf': [], 'start': [], 'end': [], 'strand': [], 'partial': []}
        orfs_aa[i] = []
        orfs_dna[i] = []
        for j, gene in enumerate(genes):
            orfs_db['contig'].append(i)
            orfs_db['orf'].append('orf_' + str(j))
            orfs_db['start'].append(gene.begin)
            orfs_db['end'].append(gene.end)
            orfs_db['strand'].append(gene.strand)
            if gene.partial_begin or gene.partial_end is True:
                orfs_db['partial'].append('yes')
            else:
                orfs_db['partial'].append('no')
            aa_record = SeqIO.SeqRecord(Seq(gene.translate()), id = i + '_orf_' + str(j), name = '', description = '')
            orfs_aa[i + '_orf_' + str(j)] = aa_record
            if gene.strand > 0:
                dna_record = SeqIO.SeqRecord(Seq(str(seq_dict[i].seq[int(gene.begin) - 1: int(gene.end)])), id = i + '_orf_' + str(j), name = '', description = '')
            else:
                dna_record = SeqIO.SeqRecord(Seq(str(seq_dict[i].seq[int(gene.begin) - 1: int(gene.end)].reverse_complement())), id = i + '_orf_' + str(j), name = '', description = '')
            orfs_dna[i + '_orf_' + str(j)] = dna_record
        orf_df = pd.DataFrame(orfs_db)
        orfs_df_list.append(orf_df)
    if len(orfs_df_list) > 1:
        orf_df = pd.concat(orfs_df_list)
    else:
        orf_df = orfs_df_list[0]
    orf_df = orf_df.reset_index(drop = True)
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return orf_df, orfs_aa, orfs_dna
