"""
functions regarding orf prediction and mapping
"""


import collections
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

import pandas as pd
import pyhmmer
import pyrodigal
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from loguru import logger
from pyhmmer.easel import Alphabet, SequenceFile
from pyhmmer.plan7 import HMM, HMMFile

from hafeZ.utils.exit import exit_error_gracefully


def get_roi_sequences(
    roi_df: pd.DataFrame, filtered_seq_dict: Dict[str, Seq], output: Path
) -> List[SeqIO.SeqRecord]:
    """
    Retrieve sequences for regions of interest (ROIs) based on the provided DataFrame and sequence data.

    Args:
        roi_df (pd.DataFrame): The DataFrame containing ROI information.
        seq (Dict[str, Seq]): A dictionary of sequences with contig IDs as keys and sequences as values.
        output (str): The folder where the ROI sequences will be written as a FASTA file.

    Returns:
        List[SeqIO.SeqRecord]: List of SeqRecord objects containing the ROI sequences.
    """

    refined_roi_fasta: Path = Path(output) / "temp_refined_roi.fasta"

    roi_list = []
    for index, row in roi_df.iterrows():
        if (row["circular"] is False) or (
            (row["circular"] is True) and (row["roi"].split("_")[-1] == "c1")
        ):
            start = int(row["start_pos"])
            end = int(row["end_pos"])
            roi_str = filtered_seq_dict[row["roi"].split("~")[0]][start : end - 1]
            roi = SeqIO.SeqRecord(
                Seq(roi_str),
                id=row["roi"],
                name="",
                description="",
            )
            roi_list.append(roi)
        # need to test the circular ones
        elif (row["circular"] is True) and (row["roi"].split("_")[-1] == "c2"):
            start_1 = int(row["start_pos"])
            end_1 = int(row["contig_len"])
            roi_1 = filtered_seq_dict[row["roi"].split("~")[0]][start_1 : end_1 - 1]
            start_2 = 0
            end_2 = int(row["end_pos"])
            roi_2 = filtered_seq_dict[row["roi"].split("~")[0]][start_2 : end_2 - 1]
            roi = SeqIO.SeqRecord(
                Seq("".join([str(seq_rec) for seq_rec in [roi_1.seq, roi_2.seq]])),
                id=row["roi"],
                name="",
                description="",
            )
            roi_list.append(roi)
    # write the ROIs
    SeqIO.write(roi_list, refined_roi_fasta, "fasta")
    return roi_list


from typing import Dict, List, Tuple

import pandas as pd
import pyrodigal
from Bio import SeqIO
from Bio.Seq import Seq


def get_orfs(
    filtered_seq_dict: Dict[str, Seq], multicontig: bool
) -> Tuple[pd.DataFrame, Dict[str, SeqIO.SeqRecord], Dict[str, SeqIO.SeqRecord]]:
    """
    Find and extract ORFs (Open Reading Frames) from a dictionary of sequences.

    Args:
        filtered_seq_dict (Dict[str, Seq]): A dictionary of sequences where keys are contig IDs and values are Bio.Seq objects.
        multicontig (bool): Indicates whether the sequences are from multiple contigs.

    Returns:
        Tuple[pd.DataFrame, Dict[str, SeqIO.SeqRecord], Dict[str, SeqIO.SeqRecord]]:
            - A pandas DataFrame containing ORF information, including contig, ORF ID, start position, end position, strand, and partial status.
            - A dictionary where keys are ORF IDs, and values are Bio.SeqRecord objects representing the corresponding amino acid sequences.
            - A dictionary where keys are ORF IDs, and values are Bio.SeqRecord objects representing the corresponding DNA sequences.
    """
    orfs_aa = {}
    orfs_dna = {}
    orfs_df_list = []
    p = pyrodigal.GeneFinder(meta=multicontig)

    for i in filtered_seq_dict:
        if multicontig is False:  # can add a length check here
            p.train(str(filtered_seq_dict[i]))
        genes = p.find_genes(str(filtered_seq_dict[i]))
        orfs_db = {
            "contig": [],
            "orf": [],
            "start": [],
            "end": [],
            "strand": [],
            "partial": [],
        }
        orfs_aa[i] = []
        orfs_dna[i] = []

        for j, gene in enumerate(genes):
            orfs_db["contig"].append(i)
            orfs_db["orf"].append("orf_" + str(j))
            orfs_db["start"].append(gene.begin)
            orfs_db["end"].append(gene.end)
            orfs_db["strand"].append(gene.strand)
            if gene.partial_begin or gene.partial_end is True:
                orfs_db["partial"].append("yes")
            else:
                orfs_db["partial"].append("no")

            aa_record = SeqIO.SeqRecord(
                Seq(gene.translate()), id=i + "_orf_" + str(j), name="", description=""
            )
            orfs_aa[i + "_orf_" + str(j)] = aa_record

            if gene.strand > 0:  # positive
                dna_record = SeqIO.SeqRecord(
                    Seq(str(filtered_seq_dict[i][int(gene.begin) - 1 : int(gene.end)])),
                    id=i + "_orf_" + str(j),
                    name="",
                    description="",
                )
            else:
                dna_record = SeqIO.SeqRecord(
                    Seq(
                        filtered_seq_dict[i][int(gene.begin) - 1 : int(gene.end)]
                    ).reverse_complement(),
                    id=i + "_orf_" + str(j),
                    name="",
                    description="",
                )
            orfs_dna[i + "_orf_" + str(j)] = dna_record

        orf_df = pd.DataFrame(orfs_db)
        orfs_df_list.append(orf_df)

    if len(orfs_df_list) > 1:
        orf_df = pd.concat(orfs_df_list)
    else:
        orf_df = orfs_df_list[0]

    orf_df = orf_df.reset_index(drop=True)

    return orf_df, orfs_aa, orfs_dna


def extract_roi_orfs(
    orf_df: pd.DataFrame,
    orf_aa: dict,
    roi_df: pd.DataFrame,
    orf_dna: dict,
    output: Path,
    min_orfs: int,
    all_zscores: bool,
    depths: Dict[Union[str, int], List[float]],
    median: float,
    mad: float,
    start_time: float,
) -> Tuple[Union[str, pd.DataFrame], Union[str, dict], Union[str, dict]]:
    """
    Extract ORF sequences for regions of interest (ROIs) based on the provided DataFrames and dictionaries.

    Args:
        orf_df (pd.DataFrame): The DataFrame containing ORF information.
        orf_aa (dict): A dictionary of ORF amino acid sequences with identifiers as keys.
        roi_df (pd.DataFrame): The DataFrame containing ROI information.
        orf_dna (dict): A dictionary of ORF DNA sequences with identifiers as keys.
        output (str): The folder where the ROI ORF sequences will be written as a FASTA file.
        min_orfs (int): The minimum number of ORFs required to keep a ROI.
        all_zscores (bool): Flag indicating whether to calculate and save z-scores.
        depths (dict): A dictionary of depth data.
        median (float): The median value.
        mad (float): The Median Absolute Deviation (MAD) value.
        start_time (float): The start time of the process.

    Returns:
        Tuple[Union[str, pd.DataFrame], Union[str, dict], Union[str, dict]]:
            - A DataFrame containing ROI information with an additional "orf_count" column.
            - A dictionary of ROI amino acid sequences with ROIs as keys.
            - A dictionary of ROI DNA sequences with ROIs as keys.
    """

    roi_aa = {}
    roi_dna = {}
    all_roi_orfs = []
    orf_no = []

    temp_aa_fasta: Path = Path(output) / "temp_aa.fasta"

    for index, row in roi_df.iterrows():
        contig = row["roi"].split("~")[0]
        counter = 1

        if (row["circular"] == False) or (
            (row["circular"] == True) and (row["roi"].split("_")[-1] == "c1")
        ):
            orfs = orf_df[
                (orf_df["contig"] == contig)
                & (orf_df["start"] >= row["start_pos"])
                & (orf_df["end"] <= row["end_pos"])
            ].copy()
        elif (row["circular"] == True) and (row["roi"].split("_")[-1] == "c2"):
            orfs_1 = orf_df[
                (orf_df["contig"] == contig)
                & (orf_df["start"] >= row["start_pos"])
                & (orf_df["end"] <= row["contig_len"])
            ].copy()
            orfs_2 = orf_df[
                (orf_df["contig"] == contig)
                & (orf_df["start"] >= 0)
                & (orf_df["end"] <= row["end_pos"])
            ].copy()
            orfs = pd.concat([orfs_1, orfs_2])

        orfs["new_orf"] = ""
        if len(orfs) < min_orfs:
            roi_df = roi_df[roi_df["roi"] != row["roi"]].copy()
        else:
            for i, r in orfs.iterrows():
                orfs.loc[i, "new_orf"] = "orf_" + str(counter)
                counter = counter + 1
            orf_no.append(len(orfs))
            roi_aa[row["roi"]] = []
            roi_dna[row["roi"]] = []
            for i in list(orfs["orf"]):
                aa = orf_aa[contig + "_" + i]
                aa.id = row["roi"] + "~" + orfs["new_orf"][orfs["orf"] == i].iloc[0]
                roi_aa[row["roi"]].append(aa)
                dna = orf_dna[contig + "_" + i]
                dna.id = row["roi"] + "~" + orfs["new_orf"][orfs["orf"] == i].iloc[0]
                roi_dna[row["roi"]].append(dna)
                all_roi_orfs.append(aa)
    if len(roi_df) < 1:
        exit_error_gracefully(output, all_zscores, depths, median, mad, start_time)
        logger.error(f"After filtering, 0 ROIs remain with at least {min_orfs} orfs. ")
    else:
        roi_df["orf_count"] = orf_no.copy()
        roi_df = roi_df.reset_index(drop=True)
        SeqIO.write(all_roi_orfs, temp_aa_fasta, "fasta")

        return roi_df, roi_aa, roi_dna


def run_pyhmmer(
    database: Path, output: Path, threads: int, evalue: float
) -> Dict[str, Tuple[str, str, float, float]]:
    """
    Runs PyHMMER to search for protein hits in a database.

    Args:
        database (Path): Path to the database directory.
        output (Path): Output directory.
        threads (int): Number of threads to use for processing.
        evalue (float): E-value threshold for PyHMMER.

    Returns:
        Dict[str, Tuple[str, str, float, float]]: A dictionary of top HMM hits for each protein,
            where keys are protein IDs, and values are tuples containing protein ID, hit ID, bit score, and E-value.
    """

    # Define the amino acid FASTA file
    amino_acid_fasta_file: Path = Path(output) / "temp_aa.fasta"

    # Define the result data structure
    Result = collections.namedtuple("Result", ["protein", "hit", "bitscore", "evalue"])

    # Run hmmscan and get all results
    results = []
    with pyhmmer.plan7.HMMFile(os.path.join(database, "all_phrogs.h3m")) as hmms:
        with pyhmmer.easel.SequenceFile(amino_acid_fasta_file, digital=True) as seqs:
            for hits in pyhmmer.hmmer.hmmscan(seqs, hmms, cpus=int(threads), E=evalue):
                protein = hits.query_name.decode()
                for hit in hits:
                    if hit.included:
                        results.append(
                            Result(protein, hit.name.decode(), hit.score, hit.evalue)
                        )

    # Get the best results for each protein
    best_results = {}
    keep_protein = set()
    for result in results:
        if result.protein in best_results:
            previous_bitscore = best_results[result.protein].bitscore
            if result.bitscore > previous_bitscore:
                best_results[result.protein] = result
                keep_protein.add(result.protein)
            elif result.bitscore == previous_bitscore:
                if best_results[result.protein].hit != result.hit:
                    keep_protein.remove(result.protein)
        else:
            best_results[result.protein] = result
            keep_protein.add(result.protein)

    return best_results


def calc_phrogs_frac(
    best_phrog_results_dict: Dict[str, Union[str, float]],
    roi_df: pd.DataFrame,
    phrog_fract: float,
    evalue: float,
    output: Path,
    all_zscores: bool,
    depths: Dict[Union[str, int], List[float]],
    median: float,
    mad: float,
    start_time: float,
) -> pd.DataFrame:
    """
    Calculate PHROGs fraction for each ROI and filter based on the threshold.

    Args:
        best_phrog_results_dict (Dict[str, Union[str, float]]): A dictionary containing the best PHROGs results
            with protein names as keys and a tuple of values including PHROG name, bitscore, and evalue as values.
        roi_df (pd.DataFrame): DataFrame containing ROI information.
        phrog_fract (float): The threshold fraction for PHROGs annotation.
        evalue (float): The E-value threshold for PHROGs annotation.
        output (Path): The output directory.
        all_zscores (bool): Flag indicating whether to calculate and save z-scores.
        depths (dict): A dictionary of depth data.
        median (float): The median value.
        mad (float): The Median Absolute Deviation (MAD) value.
        start_time (float): The start time of the process.

    Returns:
        pd.DataFrame: A filtered DataFrame containing ROI information based on the PHROGs fraction threshold.
    """

    # Function to count matching proteins with the right name (strip off orf) with evalue < evalue
    def count_matching_entries(roi_key):
        return len(
            [
                k
                for k, v in best_phrog_results_dict.items()
                if re.split(r"~orf_\d+", k)[0] == roi_key and v[3] < evalue
            ]
        )

    # Create a new column 'pyhmmer_phrog_hits' in roi_df to store the counts
    roi_df["pyhmmer_phrog_hits"] = roi_df["roi"].apply(count_matching_entries)

    # Calculate the fraction of PHROGs hits
    roi_df["frac_phrog"] = roi_df["pyhmmer_phrog_hits"] / roi_df["orf_count"]

    # Filter based on the PHROGs fraction threshold
    roi_df = roi_df[roi_df["frac_phrog"] >= phrog_fract]

    if roi_df is not None:
        if len(roi_df) > 0:
            logger.info(
                f"{roi_df.shape[0]} ROIs passed as they had more than {phrog_fract} proportion of genes annotated by PHROGs."
            )
    else:  # error 0 ROIs
        exit_error_gracefully(output, all_zscores, depths, median, mad, start_time)
        logger.error(
            f"0 ROIs passed as they had more than {phrog_fract} proportion of genes annotated by PHROGs."
        )

    return roi_df
