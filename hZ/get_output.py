import time
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import os, glob
import numpy as np
import matplotlib

def get_names(roi_df):
    print('\n{:#^50}'.format(' Finalising roi names '))
    start_time = time.time()
    counter = 1
    for index, row in roi_df.iterrows():
        roi_df.loc[index, 'roi_new'] = row['roi'].split('~')[0] + '~roi_' + str(counter)
        counter = counter + 1
    end_time = '{:.3f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return roi_df

def output_roi_seqs(roi_df,roi_dna,output_folder):
    print('\n{:#^50}'.format(' Outputting roi sequences '))
    start_time = time.time()
    final_roi_dna = []
    for i in roi_df['roi'].unique():
        for j in roi_dna:
            if j.id == i:
                j.id = roi_df['roi_new'][roi_df['roi'] == i].iloc[0]
                j.description = ''
                j.name = ''
                final_roi_dna.append(j)
    SeqIO.write(final_roi_dna, output_folder + '/hafeZ_all_roi_seqs.fasta', 'fasta')
    end_time = '{:.3f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))

def output_roi_orfs(roi_orf_dna,roi_df,output_folder,roi_orf_aa):
    print('\n{:#^50}'.format(' Outputting roi orfs '))
    start_time = time.time()
    for i in roi_df['roi'].unique():
        orf_list = []
        for j in roi_orf_dna[i]:
            name = roi_df['roi_new'][roi_df['roi'] == i].iloc[0]
            j.id = name + '~' + j.id.split('~')[-1]
            orf_list.append(j)
        SeqIO.write(orf_list, output_folder + '/hafeZ_orfs_dna_' + name + '.fasta', 'fasta')
    for i in roi_df['roi'].unique():
        orf_list = []
        for j in roi_orf_aa[i]:
            name = roi_df['roi_new'][roi_df['roi'] == i].iloc[0]
            j.id = name + '~' + j.id.split('~')[-1]
            orf_list.append(j)
        SeqIO.write(orf_list, output_folder + '/hafeZ_orfs_aa_' + name + '.faa', 'fasta')
    end_time = '{:.3f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))

def output_prophage_graphs(roi_df,depths,output_folder,median,mad):
    print('\n{:#^50}'.format(' Making roi figures '))
    start_time = time.time()
    matplotlib.use('Agg')
    roi_df['contig'] = roi_df['roi'].str.split('~').str[0]
    roi_df = roi_df.astype({'start_pos': int, 'end_pos': int})
    for i in roi_df['contig'].unique():
        z = [((0.6745*(x - median))/mad) for x in depths[i]]
        fig = plt.gcf()
        fig.set_size_inches(18.5, 10.5)
        plt.plot(z)
        pos = 4
        for index, row in roi_df[roi_df['contig'] == i].iterrows():
            if (row['circular'] == False) | ((row['circular'] == True) & (row['roi'].split('_')[-1] == 'c1')):
                plt.axvspan(row['start_pos'], row['end_pos'], color='r', alpha=0.5, lw=0)
                plt.annotate(row['roi_new'], xy=(row['start_pos'], pos), xytext=(row['start_pos'] + len(z)/10, pos + 1),
                        arrowprops=dict(facecolor='black', arrowstyle = '->', connectionstyle='arc3', lw =2))
            elif (row['circular'] == True) & (row['roi'].split('_')[-1] == 'c2'):
                plt.axvspan(row['start_pos'], row['contig_len'], color='r', alpha=0.5, lw=0)
                plt.axvspan(0,row['end_pos'], color='r', alpha=0.5, lw=0)
                plt.annotate(row['roi_new'], xy=(row['start_pos'], pos), xytext=(row['start_pos'] + len(z)/10, pos + 1),
                        arrowprops=dict(facecolor='black', arrowstyle = '->', connectionstyle='arc3', lw =2))
                plt.annotate('', xy=(row['end_pos'], pos),
                        arrowprops=dict(facecolor='black', arrowstyle = '->', connectionstyle='arc3', lw =2))
            pos = pos + (np.max(z)/len(roi_df[roi_df['contig'] == i]))
        plt.savefig(output_folder + '/hafeZ_prophage_for_' + i + '.png', format = 'png')
        plt.clf()
    end_time = '{:.3f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))

def output_hmm_table(roi_df,output_folder):
    print('\n{:#^50}'.format(' Outputting hmm hit table '))
    start_time = time.time()
    df = pd.read_csv(output_folder + '/temp_hmms.tab', sep='\t')
    df = df[df['orf'].str.split('~').str[:-1].str.join('~').isin(list(roi_df['roi']))]
    for index, row in df.iterrows():
        old = '~'.join(row['orf'].split('~')[:-1])
        df.loc[index, 'orf'] = roi_df['roi_new'][roi_df['roi'] == old].iloc[0]
    df.columns = ['vog_no','orf_name','e_value','orf_no','vog_phage_taxonomy','vog_description']
    df = df[['orf_name','orf_no','vog_no','e_value','vog_phage_taxonomy','vog_description']]
    df = df.reset_index(drop=True)
    sort = (df.assign(orf_no2=df['orf_no'].str.extract(r'(\d+)$').astype(int)).sort_values(['orf_name', 'orf_no2']).index)
    df = df.iloc[sort]
    df.to_csv(output_folder + '/hafeZ_hmm_hits.tsv', sep='\t', index=False)
    end_time = '{:.3f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))

def output_roi_table(roi_df,output_folder):
    print('\n{:#^50}'.format(' Outputting roi table '))
    start_time = time.time()
    for index, row in roi_df.iterrows():
        if row['circular'] == False:
            roi_df.loc[index, 'roi_length'] = row['end_pos'] - row['start_pos'] + 1
        elif row['circular'] == True:
            len_1 = row['contig_len'] - row['start_pos']
            len_2 = row['end_pos'] - 0
            roi_df.loc[index, 'roi_length'] = len_1 + len_2
    roi_df = roi_df[['roi_new', 'start_pos', 'start_count', 'end_pos', 'end_count', 'roi_length', 'orf_count','frac_pvog','circular','med_z','attL_seq','attR_seq', 'contig_split','longest_below_z']].copy()
    roi_df.columns = ['roi', 'start_pos', 'start_count', 'end_pos', 'end_count', 'roi_length', 'orf_count','frac_pvog','circular','med_z','attL_seq','attR_seq', 'contig_split', 'longest_below_z']
    for index, row in roi_df.iterrows():
        roi_df.loc[index, 'start_pos'] = row['start_pos'] + 1
        roi_df.loc[index, 'frac_pvog'] = round(row['frac_pvog'],2)
        roi_df.loc[index, 'roi_length'] = int(row['roi_length'])
    roi_df.to_csv(output_folder + '/hafeZ_summary_all_rois.tsv', sep='\t', index=False)
    end_time = '{:.3f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))

def output_no_roi(output_folder):
    roi_df = pd.DataFrame(data = {'roi': [np.nan], 'start_pos': [np.nan], 'end_pos': [np.nan], 'roi_length': [np.nan], 'orf_count': [np.nan], 'frac_pvog': [np.nan], 'circular': [np.nan]})
    roi_df.to_csv(output_folder + '/hafeZ_summary_all_rois.tsv', sep='\t', index=False)

def output_contig_Z(depths,output_folder,median,mad):
    print('\n{:#^50}'.format(' Making figures of contig Z-scores'))
    start_time = time.time()
    matplotlib.use('Agg')
    for i in depths:
        z = [((0.6745*(x - median))/mad) for x in depths[i]]
        fig = plt.gcf()
        fig.set_size_inches(18.5, 10.5)
        plt.plot(z)
        pos = 4
        plt.savefig(output_folder + '/zscores_for_contig' + i + '.png', format = 'png')
        plt.clf()
    end_time = '{:.3f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))

def clean_up(output_folder):
    print('\n{:#^50}'.format(' Doing final tidy '))
    start_time = time.time()
    for f in glob.glob(output_folder + '/temp*'):
        os.remove(f)
    end_time = '{:.3f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))


def output_hmm_table_phrogs(roi_df,output_folder,db_path):
    print('\n{:#^50}'.format(' Outputting hmm hit table '))
    start_time = time.time()
    df = pd.read_csv(output_folder + '/temp_hmms.tab', sep='\t', names = ['orf','vog_no','e_value','orf_no'])    
    df = df[df['orf'].str.split('~').str[:-1].str.join('~').isin(list(roi_df['roi']))]
    phrogs_db = pd.read_csv(db_path + '/phrogs_table_almostfinal_plusGO_wNA_utf8.tsv', sep='\t', usecols = [0,6],  names=['phrog','category'])
    phrogs_db['phrog'] = 'phrog_' + phrogs_db['phrog'].astype(str)
    for index, row in df.iterrows():
        old = '~'.join(row['orf'].split('~')[:-1])
        df.loc[index, 'orf'] = roi_df['roi_new'][roi_df['roi'] == old].iloc[0]
        if row['e_value'] > 0.00001:
            df.loc[index, 'vog_no'] = np.nan 
            df.loc[index, 'e_value'] = np.nan
        df.loc[index, 'number'] = int(row['orf_no'].split('_')[-1])
        df.loc[index,'vog_description'] = phrogs_db['category'][phrogs_db['phrog'] == row['vog_no']].iloc[0]
    df.columns = ['roi','vog_no','e_value','orf_no','number','vog_description']
    df = df.sort_values(by=['roi','number'])
    df = df[['roi','orf_no','vog_no','e_value','vog_description']]
    df = df.reset_index(drop=True)
    df.to_csv(output_folder + '/hafeZ_hmm_hits.tsv', sep='\t', index=False)
    end_time = '{:.3f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))