import time
import pandas as pd
from scipy.signal import savgol_filter
from scipy.stats import zscore
from scipy.signal import find_peaks
import numpy as np
import sys
import matplotlib.pyplot as plt
from statistics import mean
from statistics import stdev
import multiprocessing
from functools import partial


def savgol(df, bin_size):
    df = df[1]
    raw_dict = {}
    smoothed_dict = {}
    raw_dict[df['acc'].unique()[0]] = df['depth']
    z = savgol_filter(df['depth'],bin_size,2)
    smoothed_dict[df['acc'].unique()[0]] = z
    return [raw_dict,smoothed_dict]


def smooth_depths(output_folder, bin_size, threads):
    print('\n{:#^50}'.format(' Smoothing signal '))
    start_time = time.time()
    df = pd.read_csv(output_folder + '/temp_minimap.regions.bed.gz', compression = 'gzip', sep='\t',names=('acc','start','end','depth'))
    df = df.groupby('acc')
    pool = multiprocessing.Pool(processes=threads)
    dicts = pool.map(partial(savgol, bin_size = bin_size), df)
    raw_dict = {}
    smoothed_dict = {}
    for i in dicts:
        raw_dict[list(i[0].keys())[0]] = list(i[0].values())[0]
        smoothed_dict[list(i[1].keys())[0]] = list(i[1].values())[0]
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return smoothed_dict, raw_dict

def get_ZScores(depths, multicontig):
    print('\n{:#^50}'.format(' Getting Z-scores '))
    start_time = time.time()
    x = {}
    depths_list = np.empty(shape=0)
    for i in depths:
        depths_list = np.append(depths_list,depths[i])
    median = np.median(depths_list)
    mad = np.median(np.absolute(depths_list - median))
    for i in depths:
        x[i] = [((0.6745*(x - median))/mad) for x in depths[i]]
        x[i] = np.insert(x[i],0,0,axis=0) # add to start so peaks will find if whole contig passes threshold
        x[i] = np.append(x[i],0) # same as above but for end
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    return x, median, mad

def get_ROIs(zscores,cutoff,output_folder,width):
    print('\n{:#^50}'.format(' Finding potential regions of interest '))
    start_time = time.time()
    x = {}
    df = {'accession': [], 'start': [], 'end': [], 'length': [], 'contig_len': [], 'contig': []}
    for i in zscores:
        x[i] = {'start': [], 'end': []}
        above_threshold = np.diff(zscores[i] > cutoff, prepend=False)
        rise_above_thres = np.argwhere(above_threshold)[::2,0]
        sink_below_thres = np.argwhere(above_threshold)[1::2,0]
        contig_len = len(np.array(zscores[i]))
        x[i]['start'] = rise_above_thres
        x[i]['end'] = sink_below_thres
        for j in np.arange(0,len(x[i]['start'])):
            if x[i]['end'][j]- x[i]['start'][j] > width:##### added this line
                df['length'].append(x[i]['end'][j]- x[i]['start'][j])
                df['accession'].append('{}~{}~{}~{}'.format(i,x[i]['start'][j],x[i]['end'][j],j))
                df['start'].append(x[i]['start'][j])
                df['end'].append(x[i]['end'][j])
                df['contig_len'].append(contig_len)
                df['contig'].append(i.split('~')[0])
    df = pd.DataFrame(df)
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    if len(df) < 1:
        return 'exit'
    else:
        return df
