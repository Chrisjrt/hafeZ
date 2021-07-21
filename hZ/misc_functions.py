import pandas as pd
import multiprocessing
import numpy as np

def split_df_for_multiproc(df,threads):
    chunk_size = int(df.shape[0]/threads)
    if chunk_size == 0:
        chunk_size = 1
    chunks = [df.iloc[i:i + chunk_size,:] for i in range(0, df.shape[0], chunk_size)]
    pool = multiprocessing.Pool(processes=threads)
    return chunks, pool
