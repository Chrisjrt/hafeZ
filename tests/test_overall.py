"""
Unit tests for hafeZ overall

From the hafeZ directory:
Usage: pytest .

"""

# import
import os
import shutil

# import functions
import subprocess
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import pytest
from loguru import logger

from hafeZ.utils.util import remove_directory

# import functions


# test data
test_data = Path("tests/test_data")
short_data = Path(f"{test_data}/short")
long_data = Path(f"{test_data}/long")
database_dir = Path(f"{test_data}/database")
logger.add(lambda _: sys.exit(1), level="ERROR")
threads = 4


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


temp_dir = Path(f"{test_data}/fake_out")


def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """executes shell command and returns stdout if completes exit code 0
    Parameters
    ----------
    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""

    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"FAILED: {cmnd}\n{err}")
    return out.decode("utf8") if out is not None else None


def test_download(tmp_dir):
    """test hafeZ database"""
    cmd = f"hafeZ database -d {database_dir}"
    exec_command(cmd)


def test_short(tmp_dir):
    """test hafeZ short"""
    input_r1: Path = f"{short_data}/combo_R1.fastq.gz"
    input_r2: Path = f"{short_data}/combo_R2.fastq.gz"
    input_genome: Path = f"{short_data}/C222.fasta"
    coverage: float = 3.5  # to test subsampling
    cmd = f"hafeZ short -g {input_genome} -1 {input_r1} -2 {input_r2} -o {tmp_dir} -t {threads} -S -C {coverage} -d {database_dir} --all_zscores -f"
    exec_command(cmd)


def test_long(tmp_dir):
    """test hafeZ long"""
    input_longreads: Path = f"{long_data}/combo.fastq.gz"
    input_genome: Path = f"{short_data}/C222.fasta"
    coverage: float = 3.5  # to test subsampling
    cmd = f"hafeZ long -g {input_genome} -l {input_longreads} -o {tmp_dir} -t {threads} -S -C {coverage} -d {database_dir} --all_zscores -f"
    exec_command(cmd)


class testFails(unittest.TestCase):
    """Tests for hafeZ fails"""

    def test_long_fail_no_longreads(self):
        with self.assertRaises(RuntimeError):
            """test hafeZ long if no longreads are input as a mock test case"""
            input_longreads: Path = f"{long_data}/combo.fastq.gz"
            input_genome: Path = f"{short_data}/C222.fasta"
            coverage: float = 3.5  # to test subsampling
            cmd = f"hafeZ long -g {input_genome} -o {temp_dir} -t {threads} -S -C {coverage} -d {database_dir} --all_zscores -f"
            exec_command(cmd)


remove_directory(temp_dir)
remove_directory(f"{database_dir}")
