import urllib.request
import sys
import time
import os
import tarfile
import shutil
import pandas as pd
import glob
import subprocess
import numpy as np

def progress(count, block_size, total_size): # code for progress bar taken from https://blog.shichao.io/2012/10/04/progress_speed_indicator_for_urlretrieve_in_python.html
    global start_time
    if count == 0:
        start_time = time.time()
        return
    duration = time.time() - start_time
    progress_size = int(count * block_size)
    speed = int(progress_size / (1024 * duration))
    percent = int(count * block_size * 100 / total_size)
    sys.stdout.write("\rProgress: %d%%, Downloaded: %d MB, Speed: %d KB/s" %
                    (percent, progress_size / (1024 * 1024), speed))
    sys.stdout.flush()

def get_file(url,filename,db_folder):
    filename = db_folder + '/' + filename
    urllib.request.urlretrieve(url, filename, progress)

def unpack_dbs(db_folder,filename):
    tar = tarfile.open(db_folder + '/' + filename, "r:gz")
    tar.extractall(path = db_folder + '/')
    tar.close()
    os.remove(db_folder + '/' + filename)

def get_dbs(db_folder):
    print('\n{:#^50}'.format(' Downloading pVOGs files '))
    start_time = time.time()
    os.makedirs(db_folder, exist_ok=False)
    HMMs = 'http://dmk-brain.ecn.uiowa.edu/pVOGs/downloads/All/AllvogHMMprofiles.tar.gz' # url for pVOGs HMMs
    VOG_tables = 'http://dmk-brain.ecn.uiowa.edu/pVOGs/downloads/All/Allvogtables.tar.gz' # url for pVOGs tables
    urls = [HMMs, VOG_tables]
    counter = 1
    for i in urls:
        print('\nGetting file {}/{}:\n'.format(counter,len(urls)))
        get_file(i,i.split("/")[-1],db_folder)
        print('\nunpacking the file\n')
        unpack_dbs(db_folder,i.split("/")[-1])
        counter = counter + 1
    os.system("cat {}/AllvogHMMprofiles/VOG*.hmm > {}/combined_AllvogHMMprofiles.hmm".format(db_folder,db_folder))
    shutil.rmtree(db_folder + '/AllvogHMMprofiles')
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('\n{:#^50}'.format(' Done: ' + end_time + ' seconds '))

def process_vogsTable(db_folder):
    print('\n{:#^50}'.format(' Processing pVOGs table '))
    start_time = time.time()
    df_list = []
    for f in glob.glob(db_folder + '/Allvogtables/VOG*.txt'):
        with open(f,'r') as i:
            VOG = i.readline().split(':')[0][1:]
        df = pd.read_csv(f,sep='\t',comment='#', names = ('type','phage_name','genome_acc','prot_acc','coords','length','description'))
        df['vog'] = VOG
        df_list.append(df)
    df_list = pd.concat(df_list).reset_index(drop=True)
    df_list = df_list.groupby('vog')[['type','phage_name', 'description']].agg(lambda x: x.tolist()).reset_index()
    df_list['description'] = df_list.vog.map(df_list.groupby('vog')['description'].sum().apply(np.unique))
    df_list['type'] = df_list.vog.map(df_list.groupby('vog')['type'].sum().apply(np.unique))
    shutil.rmtree(db_folder + '/Allvogtables/')
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('\n{:#^50}'.format(' Done: ' + end_time + ' seconds '))
    df_list.to_csv(db_folder + '/combined_Allvogtables.txt')

def process_vogsHMM(db_folder):
    print('\n{:#^50}'.format(' Pressing pVOGs HMMs '))
    start_time = time.time()
    db_folder = db_folder + '/combined_AllvogHMMprofiles.hmm'
    subprocess.run(['hmmpress',
                    db_folder])
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('\n{:#^50}'.format(' Done: ' + end_time + ' seconds '))
