import collections
import os
import random
import string
from cmath import nan
from pathlib import Path
from re import T

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import GC
from loguru import logger
from util import remove_directory, remove_file, touch_file
from typing import Dict, Tuple, Union


class Haf:
    """hafez Class"""

    def __init__(
        self,
        output: Path = "output_hafeZ/",
        fasta: Path = "fasta",
        database: Path = "database/",
        # threads: int = 1,
        prefix: str = "hafeZ",
        filtered_seq_dict: Dict[str, str] = {"contig1": "AGT"},
        multicontig: bool = True,
        fasta_filtered: Path = "output_hafeZ/temp.fasta",
        raw_seq_dict: Dict[str, str] = {"contig1": "AGT"},
    ) -> None:
        """
        Parameters
        --------
        output: Path,
            output directory
        fasta: Path,
            genome FASTA
        database: Path,
            hafeZ database directory
        prefix: str, prefix for output
            prefix
        filtered_seq_dict: Dict[str, str]
            filtered contig dictionary
        multicontig: bool
            True if genome has >1 passing contig
        fasta_filtered: Path
            filtered genom
        raw_seq_dict:  Dict[str, str]
            raw contig dictionary
        """
        self.output = output
        self.fasta = fasta
        self.database = database
        self.prefix = prefix
        self.filtered_seq_dict = filtered_seq_dict
        self.multicontig = multicontig
        self.fasta_filtered = fasta_filtered
        self.raw_seq_dict = raw_seq_dict
