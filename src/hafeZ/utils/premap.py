import time
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from loguru import logger
from pathlib import Path
from typing import Dict, Tuple, Union
from hafeZ.utils.util import count_bases_in_fastq
from hafeZ.utils.external_tools import ExternalTool
from hafeZ.utils.exit import exit_error_gracefully_premap

class Premap:
    """Premapping Class for hafeZ"""

    def __init__(
        self,
        filtered_seq_dict: Dict[str, str] = {"contig1": "AGT"},
        multicontig: bool = True,
        fasta_filtered: Path = "output_hafeZ/temp.fasta",
        raw_seq_dict: Dict[str, str] = {"contig1": "AGT"},
        genome_length: int = 0,
        coverage: float = 100.0,
        estimated_coverage: float = 100.0,
        output: Path = "output_hafeZ",
        start_time: float = 0.0
    ) -> None:
        """
        Parameters
        --------
        filtered_seq_dict: Dict[str, str]
            filtered contig dictionary
        multicontig: bool
            True if genome has >1 passing contig
        fasta_filtered: Path
            filtered genom
        raw_seq_dict:  Dict[str, str]
            raw contig dictionary
        genome_length: int
            genome length of the fasta_filtered
        coverage: float
            given desired coverage
        estimated_coverage: float
            calculated estimated coverage
        output: Path
            output directory
        start_time: float
            start time
        """
        self.filtered_seq_dict = filtered_seq_dict
        self.multicontig = multicontig
        self.fasta_filtered = fasta_filtered
        self.raw_seq_dict = raw_seq_dict
        self.genome_length = genome_length
        self.coverage = coverage
        self.estimated_coverage = estimated_coverage
        self.output = output
        self.start_time = start_time

    def process_fasta(self, genome: Path, min_contig_len: int) -> None:
        """
        Takes an input genome FASTA file, keeps only contigs with a length greater than min_contig_len, and returns relevant data.

        Parameters:
            genome (Path): Path to the input FASTA file containing biological sequences.
            output (Path): Path to the output directory where the filtered FASTA file will be saved.
            min_contig_len (int): Minimum length threshold for retaining contigs.

        Returns:
            None

        Saves in class
        seq_dict (Dict[str, str]): A dictionary containing filtered contigs with sequence IDs as keys and sequences as values.
        multicontig (bool): A boolean indicating whether there are multiple contigs in the input FASTA file.
        seq_temp_out (Path): Path to the temporary filtered FASTA file.
        raw_seq_dict (Dict[str, str]): A dictionary containing all sequences (including those filtered out) with sequence IDs as keys and sequences as values.
        """
        logger.info(f"Processing genome FASTA {genome}")

        raw_seq_dict = {}
        filtered_seq_dict = {}
        seq_temp_out: Path = self.output / "temp_genome.fasta"
        seq_list = []

        # keeps only contigs > min_contig_len
        # Iterate over the FASTA file and process sequences
        for record in SeqIO.parse(genome, "fasta"):
            raw_seq_dict[record.id] = str(record.seq)
            if len(record.seq) > min_contig_len:
                filtered_seq_dict[record.id] = str(record.seq)
                seq_list.append(record)
        # check how many contigs passed the thresholds
        if len(filtered_seq_dict) < 1:

            exit_error_gracefully_premap(self.output, self.start_time)
            logger.error(
                f"No contigs were over {min_contig_len} that passed initial filtering. Please check your input FASTA {genome}."
            )
            
        elif len(filtered_seq_dict) == 1:
            multicontig = False
        elif len(filtered_seq_dict) > 1:
            multicontig = True
        # write to file as temp FASTA
        SeqIO.write(seq_list, seq_temp_out, "fasta")

        # save in the class
        self.filtered_seq_dict = filtered_seq_dict
        self.multicontig = multicontig
        self.seq_temp_out = seq_temp_out
        self.raw_seq_dict = raw_seq_dict

    def get_genome_length(self) -> None:
        """
        Calculate the total length of the genome by summing the lengths of all sequences in raw_seq_dict.

        Parameters:
            None
        
        Returns:
            None

        This method calculates the total length of the genome by iterating through the sequences
        in the 'raw_seq_dict' dictionary and summing their lengths. The result is stored as the
        'genome_length' attribute of the object for later use.
        """
        logger.info("Calculating genome length.")
        genome_length = 0
        for i in self.raw_seq_dict:
            genome_length += len(self.raw_seq_dict[i])
        self.genome_length = genome_length


    def calculate_estimated_coverage_short(self, reads1: Path, reads2: Path) ->  None :
        """
        Calculate and estimate the predicted overall coverage based on input read files.

        Parameters:
            reads1 (Path): Path to the first input FASTQ file containing sequencing reads.
            reads2 (Path): Path to the second input FASTQ file containing sequencing reads.

        Returns:
            None

        This method calculates the estimated overall coverage of a genome based on the provided
        input read files ('reads1' and 'reads2'). It first counts the total number of bases in
        the input FASTQ files and then divides this total by the genome length to estimate the
        coverage. The result is stored as the 'estimated_coverage' attribute of the object.
        """

        # calculate the coverage
        

        file1_bases = count_bases_in_fastq(reads1)
        file2_bases = count_bases_in_fastq(reads2)
        total_bases = file1_bases + file2_bases

        estimated_coverage = total_bases / self.genome_length
        self.estimated_coverage = estimated_coverage

    def subsample_short(self, reads1: Path, reads2: Path, logdir: Path ) ->  Tuple[Path, Path] :
        """
        Subsample sequencing reads to achieve a target coverage if necessary.

        Parameters:
            reads1 (Path): Path to the input FASTQ file containing the first set of sequencing reads.
            reads2 (Path): Path to the input FASTQ file containing the second set of sequencing reads.
            logdir (Path): Path to the directory where log files will be stored.

        Returns:
            Tuple[Path, Path]: A tuple containing the paths to the subsampled first and second read files.

        This method subsamples sequencing reads to achieve a target coverage if the estimated coverage
        ('self.estimated_coverage') is greater than the desired coverage ('self.coverage'). It uses the
        'rasusa' tool to perform subsampling. If the estimated coverage is less than or equal to the desired
        coverage, no subsampling is performed, and the original read paths are returned.

        If subsampling is required, the subsampled read files are stored in the 'output' directory with names
        "subsampled_R1.fastq.gz" and "subsampled_R2.fastq.gz" for the first and second read sets, respectively.

        Note: This method assumes that the 'rasusa' tool is available and properly configured.

        """
        
        # set path of reads to be used later
        # will be unchanged if no subsample
        outreads1: Path = self.output / "subsampled_R1.fastq.gz"
        outreads2: Path = self.output / "subsampled_R2.fastq.gz"

        # runs rasusa if sampling
        if self.estimated_coverage > self.coverage:
            logger.info(f"Estimated coverage {self.estimated_coverage:.2f}x is more than desired coverage ({self.coverage}x).")
            logger.info(f"Subsampling with rasusa.")
            rasusa = ExternalTool(
            tool="rasusa",
            input=f"",
            output=f"",
            params=f'-i {reads1} -i {reads2} --coverage {str(self.coverage)} -s 13  --genome-size {str(self.genome_length)} -o {outreads1} -o {outreads2} -O G',
            logdir=logdir,
        )
            
            ExternalTool.run_tool(rasusa)

        elif self.coverage > self.estimated_coverage:
            logger.info(f"Estimated coverage {self.estimated_coverage:.2f} is less than desired coverage ({self.coverage}x).")
            logger.info(f"Continuing with estimated coverage {self.estimated_coverage:.2f}")
            outreads1 = reads1
            outreads2 = reads2

        return outreads1, outreads2


             
             

        

