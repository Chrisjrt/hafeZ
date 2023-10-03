"""MISC FUNCTIONS
"""

import os
import binascii
import multiprocessing
import shutil
import gzip
import sys
import time
import click
from loguru import logger


class OrderedCommands(click.Group):
    """This class will preserve the order of subcommands, which is useful when printing --help"""

    def list_commands(self, ctx: click.Context):
        return list(self.commands)


# function to touch create a file
# https://stackoverflow.com/questions/12654772/create-empty-file-using-python
def touch_file(path):
    with open(path, "a"):
        os.utime(path, None)


def remove_file(file_path):
    if os.path.exists(file_path):
        os.remove(file_path)


def remove_directory(dir_path):
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)


def hafeZ_base(rel_path):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def get_version():
    with open(hafeZ_base("VERSION"), "r") as f:
        version = f.readline()
    return version


def echo_click(msg, log=None):
    click.echo(msg, nl=False, err=True)
    if log:
        with open(log, "a") as lo:
            lo.write(msg)


def print_citation():
    with open(hafeZ_base("CITATION"), "r") as f:
        for line in f:
            echo_click(line)


log_fmt = (
    "[<green>{time:YYYY-MM-DD HH:mm:ss}</green>] <level>{level: <8}</level> | "
    "<level>{message}</level>"
)

"""
logo
"""


def print_splash():
    click.echo(
        """\b

    /\  /\ 
 _  \/  \/    __      _____
| |__   __ _ / _| ___|__  /
| '_ \ / _` | |_ / _ \ / / 
| | | | (_| |  _|  __// /___ 
|_| |_|\__,_| |  \___/_____/
            | |
            |/
"""
    )


"""
begin function
"""


def begin_hafeZ(params):
    """
    begins hafeZ
    params: params a dictionary of params for hazed
    returns start time
    """
    # get start time
    start_time = time.time()
    # initial logging stuff
    log_file = os.path.join(params["--output"], f"hafeZ_{start_time}.log")
    # adds log file
    logger.add(log_file)
    logger.add(lambda _: sys.exit(1), level="ERROR")

    print_splash()

    logger.info("hafeZ: Identify inducible prophages through bacterial genomic read mapping.")
    logger.info(f"You are using hafeZ version {get_version()}")
    logger.info(f"Repository homepage is https://github.com/Chrisjrt/hafeZ.")
    logger.info(f"Listing parameters.")

    for key, value in params.items():
        logger.info(f"Parameter: {key} {value}.")
    return start_time





"""
from phiSpy
https://github.com/linsalrob/PhiSpy/blob/2e0e96ad3c5b92a2154b18de37d4c2865687c06e/PhiSpyModules/helper_functions.py#L46
"""

def is_gzip_file(f):
    """
    This is an elegant solution to test whether a file is gzipped by reading the first two characters.
    I also use a version of this in fastq_pair if you want a C version :)
    See https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed for inspiration
    :param f: the file to test
    :return: True if the file is gzip compressed else false
    """
    with open(f, 'rb') as i:
        return binascii.hexlify(i.read(2)) == b'1f8b'
    



def count_bases_in_fastq(fastq_file):
    total_bases = 0

    # if gzip
    if is_gzip_file(fastq_file) is True:
        with gzip.open(fastq_file, 'rt') as file:
            for line_number, line in enumerate(file, start=1):
                if line_number % 4 == 2:
                    total_bases += len(line.strip())
    else:
        with open(fastq_file, 'r') as file:
            for line_number, line in enumerate(file, start=1):
                if line_number % 4 == 2:
                    total_bases += len(line.strip())

    return total_bases


def split_df_for_multiproc(df, threads):
    chunk_size = int(df.shape[0] / threads)
    if chunk_size == 0:
        chunk_size = 1
    chunks = [df.iloc[i : i + chunk_size, :] for i in range(0, df.shape[0], chunk_size)]
    pool = multiprocessing.Pool(processes=threads)
    return chunks, pool
