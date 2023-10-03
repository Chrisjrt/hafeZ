"""
lots of this code taken from https://github.com/gbouras13/dnaapler
"""

import re
import shutil
import subprocess as sp
import sys
from pathlib import Path
import pyrodigal

from Bio import SeqIO
from loguru import logger


def instantiate_dirs(output_dir: str, force: bool):
    """Checks the output directory
    :param out_dir: output directory path
    :param force: force flag
    :param logger: logger
    :return: out_dir: final output directory
    """

    # Checks the output directory
    # remove outdir on force
    logger.add(lambda _: sys.exit(1), level="ERROR")
    logger.info(f"Checking the output directory {output_dir}.")
    if force is True:
        if Path(output_dir).exists():
            shutil.rmtree(output_dir)
        else:
            logger.info(
                "--force was specified even though the output directory does not already exist. Continuing."
            )
    else:
        if Path(output_dir).exists():
            logger.error(
                "Output directory already exists and force was not specified. Please specify -f or --force to overwrite the output directory."
            )

    # instantiate outdir
    if Path(output_dir).exists() is False:
        Path(output_dir).mkdir(parents=True, exist_ok=True)


def validate_fasta(input_fasta: Path):
    """
    Validates  FASTA input - checks that the input is a FASTA
    """
    logger.info(
        f"Checking that the input file {input_fasta} is in FASTA format."
    )
    # to get extension
    with open(input_fasta, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        if any(fasta):
            logger.info(f"{input_fasta} file checked.")
        else:
            logger.error(
                f"Error: {input_fasta} file is not in the FASTA format. Please check your input file."
            )


def check_dependencies():
    """Checks the versions of all dependencies
    :return:
    """

    # samtools
    try:
        process = sp.Popen(["samtools", "--version"], stdout=sp.PIPE, stderr=sp.PIPE)
        samtools_out, _ = process.communicate()
        samtools_out = samtools_out.decode()
        samtools_version = samtools_out.split("\n")[0].split(" ")[
            1
        ]  # get second line, and then second component of line
        message = f"Samtools v{samtools_version} found."
        logger.info(message)
    except Exception:
        logger.error("Samtools not found.")

    # mosdepth
    try:
        process = sp.Popen(["mosdepth", "--version"], stdout=sp.PIPE, stderr=sp.PIPE)
        samtools_out, _ = process.communicate()
        samtools_out = samtools_out.decode()
        samtools_version = samtools_out.split("\n")[0].split(" ")[
            1
        ]  # get second line, and then second component of line
        message = f"mosdepth v{samtools_version} found."
        logger.info(message)
    except Exception:
        logger.error("mosdepth not found.")

    # minimap2
    try:
        process = sp.Popen(["minimap2", "--version"], stdout=sp.PIPE, stderr=sp.PIPE)
        minimap2_out, _ = process.communicate()
        minimap2_version = minimap2_out.decode()
        minimap2_version = minimap2_version.split("\n")[0]
        message = f"minimap2 v{minimap2_version} found."
        logger.info(message)
    except Exception:
        logger.error("minimap2 not found.")

    #######
    # pyrodigal
    #######

    pyrodigal_version = pyrodigal.__version__
    pyrodigal_major_version = int(pyrodigal_version.split(".")[0])

    if pyrodigal_major_version < 3:
        logger.error("Pyrodigal is the wrong version. Please re-install pharokka.")

    logger.info(f"Pyrodigal version is v{pyrodigal_version}")
    logger.info(f"Pyrodigal version is ok.")

    # all dependencies found
    logger.info("All dependencies found.")


def is_scientific_notation(evalue):
    """
    checks if evalue is scientific notation
    """
    # Define the regular expression pattern for scientific notation
    scientific_pattern = r"^[+\-]?(\d+(\.\d*)?|\.\d+)([eE][+\-]?\d+)?$"

    # Check if the number matches the scientific notation pattern
    return bool(re.match(scientific_pattern, evalue))


def is_numeric(evalue):
    """
    checks if evalue is numeric
    """
    try:
        float(evalue)  # Attempt to convert the value to a float
        return True
    except ValueError:
        return False


def check_memory_limit(memory_limit):
    """
    checks if memory limit is a string that ends with K, G or M.
    """

    memory_limit = memory_limit.upper()
    # -1 gets the last char
    if not any(x in memory_limit[-1] for x in ["K", "G", "M"]):
        logger.error(
            f"Memory limit format must be in format of e.g. 2GB or 2MB (GB being gigabytes of RAM and MB being megabytes of RAM)"
        )
