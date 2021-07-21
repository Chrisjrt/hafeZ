import subprocess
import os
import time

def get_rough_coverage_prediction(reads1, reads2, genome_length):
    print('\n{:#^50}'.format(' Estimating predicted overall coverage '))
    start_time = time.time()
    read_bases = subprocess.run(['reformat.sh',
                        'in1=' + reads1,
                        'in2=' + reads2], capture_output=True)
    read_bases = int(read_bases.stderr.decode("utf-8").split('\t')[-4].split()[0])
    coverage = read_bases/genome_length
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + str(end_time) + ' seconds '))
    return coverage


def take_read_sample(reads1, reads2, coverage, output_folder, desired_coverage):
    print('\n{:#^50}'.format(' Subsampling reads to get desired coverage '))
    start_time = time.time()
    proportion = desired_coverage/coverage
    read1_path = output_folder + '/temp_read1.fastq.gz'
    read2_path = output_folder + '/temp_read2.fastq.gz'
    subprocess.run(['reformat.sh',
                        'in1=' + reads1,
                        'in2=' + reads2,
                        'out1=' + read1_path,
                        'out2=' + read2_path,
                        'samplerate=' + str(proportion),
                        'sampleseed=1',
                        'overwrite=true'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + str(end_time) + ' seconds '))
    return read1_path, read2_path

def minimap_paired(reads1, reads2, fasta, output_folder, threads):
    print('\n{:#^50}'.format(' Mapping reads to assembly '))
    start_time = time.time()
    db = output_folder + '/' + os.path.basename(os.path.splitext(fasta)[0])
    subprocess.run(['minimap2',
                        '-ax',
                        'sr',
                        fasta,
                        reads1,
                        reads2,
                        '-t',
                        str(threads),
                        '-o',
                        output_folder + '/temp_minimap.sam'])
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + str(end_time) + ' seconds '))

def get_bam(output_folder,threads,memory):
    print('\n{:#^50}'.format(' Generating BAM file '))
    start_time = time.time()
    print('Running samtools')
    subprocess.run(['samtools',
                        'view',
                        '-S',
                        '-b',
                        '-@',
                        str(threads),
                        '-o',
                        output_folder + '/temp_minimap.bam',
                        output_folder + '/temp_minimap.sam'])
    subprocess.run(['samtools',
                        'sort',
                        '-@',
                        str(threads),
                        '-m',
                        memory,
                        '-o',
                        output_folder + '/temp_minimap_sorted.bam',
                        output_folder + '/temp_minimap.sam'])
    print('finished running samtools')
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + str(end_time) + ' seconds '))

def get_cov(output_folder,threads):
    print('\n{:#^50}'.format(' Generating coverage depths '))
    start_time = time.time()
    print('Running mosdepth')
    subprocess.run(['samtools',
                        'index',
                        '-@',
                        str(threads),
                        '-b',
                        output_folder + '/temp_minimap_sorted.bam'])
    subprocess.run(['mosdepth',
                        '-t',
                        str(threads),
                        '-b',
                        '1',
                        '-n',
                        output_folder + '/temp_minimap',
                        output_folder + '/temp_minimap_sorted.bam'])
    print('finished running mosdepth')
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + str(end_time) + ' seconds '))

def map_clipped_reads(output_folder,threads,fasta,memory_limit):
    print('\n{:#^50}'.format(' Remapping clipped ends '))
    start_time = time.time()
    subprocess.run(['minimap2',
                    '-ax',
                    'sr',
                    fasta,
                    output_folder + '/temp_soft_clipped.fastq',
                    '-t',
                    str(threads),
                    '-o',
                    output_folder + '/temp_minimap_clipped.sam'])
    subprocess.run(['samtools',
                    'sort',
                    '-@',
                    str(threads),
                    '-m',
                    memory_limit,
                    '-o',
                    output_folder + '/temp_minimap_clipped_sorted.bam',
                    output_folder + '/temp_minimap_clipped.sam'])
    subprocess.run(['samtools',
                    'index',
                    '-@',
                    str(threads),
                    '-b',
                    output_folder + '/temp_minimap_clipped_sorted.bam'])
    end_time = '{:.2f}'.format(time.time() - start_time)
    print('{:#^50}'.format(' Done: ' + end_time + ' seconds '))
