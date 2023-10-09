import collections
import glob
import os
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from loguru import logger

from hafeZ.utils.external_tools import ExternalTool


def get_names(roi_df: pd.DataFrame) -> pd.DataFrame:
    """
    Generate new ROI names and update the DataFrame.

    Args:
    roi_df (pd.DataFrame): DataFrame containing ROI information.

    Returns:
    pd.DataFrame: Updated DataFrame with new ROI names.
    """

    counter = 1
    roi_df["roi_new"] = ""
    for index, row in roi_df.iterrows():
        roi_df.loc[index, "roi_new"] = row["roi"].split("~")[0] + "~roi_" + str(counter)
        counter += 1
    return roi_df


def get_att(
    roi_df: pd.DataFrame, filtered_seq_dict: dict, output: Path, logdir: Path
) -> pd.DataFrame:
    """
    Extracts attL and attR sequences for each ROI in the DataFrame.

    Args:
    roi_df (pd.DataFrame): DataFrame containing ROI information.
    filtered_seq_dict (dict): Dictionary of filtered sequences.
    output (Path): Path to the output directory.
    logdir (Path): Path to the log directory.

    Returns:
    pd.DataFrame: Updated DataFrame with attL and attR sequences.
    """

    left_fasta: Path = Path(output) / "temp_left.fasta"
    right_fasta: Path = Path(output) / "temp_right.fasta"
    blast_output: Path = Path(output) / "attr_blast.txt"

    left_list = []
    right_list = []
    roi_df["start_pos"] = roi_df["start_pos"].astype(int)
    roi_df["end_pos"] = roi_df["end_pos"].astype(int)

    # instantiate the attL and attR
    roi_df["attL_seq"] = ""
    roi_df["attR_seq"] = ""

    for index, row in roi_df.iterrows():
        if not pd.isna(row["contig_split"]):
            roi_df.loc[index, "attL_seq"] = np.nan
            roi_df.loc[index, "attR_seq"] = np.nan
        if row["circular"] is True:
            roi_df.loc[index, "attL_seq"] = np.nan
            roi_df.loc[index, "attR_seq"] = np.nan
        else:
            left = filtered_seq_dict[row["contig"]][
                row["start_pos"] - 100 : row["start_pos"] + 100
            ]
            right = filtered_seq_dict[row["contig"]][
                row["end_pos"] - 100 : row["end_pos"] + 100
            ]

            left_area = SeqIO.SeqRecord(
                Seq(left),
                id=row["roi"],
                name="",
                description="",
            )

            right_area = SeqIO.SeqRecord(
                Seq(right),
                id=row["roi"],
                name="",
                description="",
            )

            left_list.append(left_area)
            right_list.append(right_area)
    if len(left_list) > 0 and len(right_list) > 0:
        # write to file
        SeqIO.write(left_list, left_fasta, "fasta")
        SeqIO.write(right_list, right_fasta, "fasta")

        blast = ExternalTool(
            tool="blastn",
            input=f"",
            output=f"",
            params=f' -query {left_fasta} -subject {right_fasta} -out {blast_output} -evalue 1 -task blastn-short -outfmt  "6 qseqid qstart qend sseqid sstart send evalue qseq sseq length" ',
            logdir=logdir,
        )  # healthy e value of 1

        ExternalTool.run_tool(blast)

        blast_df = pd.read_csv(
            blast_output,
            sep="\t",
            names=(
                "qname",
                "qstart",
                "qend",
                "sname",
                "sstart",
                "send",
                "e",
                "qseq",
                "sseq",
                "length",
            ),
        )
        # get only hits > 11 bp
        blast_df = blast_df.astype({"e": float})
        blast_df = blast_df.astype({"qseq": str})
        blast_df = blast_df[blast_df["length"] > 11]
    else:
        column_names = [
            "qname",
            "qstart",
            "qend",
            "sname",
            "sstart",
            "send",
            "e",
            "qseq",
            "sseq",
            "length",
        ]
        blast_df = pd.DataFrame(columns=column_names)

    # add attL and attR to dataframe
    for index, row in roi_df.iterrows():
        sub_df = blast_df[
            (blast_df["qname"] == row["roi"]) & (blast_df["sname"] == row["roi"])
        ]
        if len(sub_df) > 0:
            # get lowest d value
            sub_df = sub_df.sort_values(by=["e"]).copy()
            roi_df.loc[index, "attL_seq"] = str(sub_df["qseq"].iloc[0])
            roi_df.loc[index, "attR_seq"] = str(sub_df["sseq"].iloc[0])
        else:
            roi_df.loc[index, "attL_seq"] = np.nan
            roi_df.loc[index, "attR_seq"] = np.nan
    return roi_df


def output_roi_seqs(roi_df, roi_dna, output):
    """
    Writes ROI DNA sequences to a FASTA file.

    Args:
    roi_df (pd.DataFrame): DataFrame containing ROI information.
    roi_dna (list): List of ROI DNA sequences as SeqRecords.
    output (Path): Path to the output directory.
    """
    final_roi_dna = []
    for i in roi_df["roi"].unique():
        for j in roi_dna:
            if j.id == i:
                j.id = roi_df["roi_new"][roi_df["roi"] == i].iloc[0]
                j.description = ""
                j.name = ""
                final_roi_dna.append(j)
    SeqIO.write(final_roi_dna, output / "hafeZ_all_roi_seqs.fasta", "fasta")


def output_prophage_graphs(
    roi_df: pd.DataFrame, depths: dict, output: Path, median: float, mad: float
) -> None:
    """
    Generates and saves prophage coverage graphs for each contig.

    Args:
    roi_df (pd.DataFrame): DataFrame containing ROI information.
    depths (dict): Dictionary of coverage depths for each contig.
    output_folder (Path): Path to the output directory.
    median (float): Median depth for normalization.
    mad (float): Median Absolute Deviation (MAD) for normalization.
    """
    matplotlib.use("Agg")
    roi_df["contig"] = roi_df["roi"].str.split("~").str[0]
    roi_df = roi_df.astype({"start_pos": int, "end_pos": int})
    for i in roi_df["contig"].unique():
        z = [((0.6745 * (x - median)) / mad) for x in depths[i]]
        fig = plt.gcf()
        fig.set_size_inches(18.5, 10.5)
        plt.plot(z)
        pos = 4
        for index, row in roi_df[roi_df["contig"] == i].iterrows():
            if (row["circular"] == False) | (
                (row["circular"] == True) & (row["roi"].split("_")[-1] == "c1")
            ):
                plt.axvspan(
                    row["start_pos"], row["end_pos"], color="r", alpha=0.5, lw=0
                )
                plt.annotate(
                    row["roi_new"],
                    xy=(row["start_pos"], pos),
                    xytext=(row["start_pos"] + len(z) / 10, pos + 1),
                    arrowprops=dict(
                        facecolor="black", arrowstyle="->", connectionstyle="arc3", lw=2
                    ),
                )
            elif (row["circular"] == True) & (row["roi"].split("_")[-1] == "c2"):
                plt.axvspan(
                    row["start_pos"], row["contig_len"], color="r", alpha=0.5, lw=0
                )
                plt.axvspan(0, row["end_pos"], color="r", alpha=0.5, lw=0)
                plt.annotate(
                    row["roi_new"],
                    xy=(row["start_pos"], pos),
                    xytext=(row["start_pos"] + len(z) / 10, pos + 1),
                    arrowprops=dict(
                        facecolor="black", arrowstyle="->", connectionstyle="arc3", lw=2
                    ),
                )
                plt.annotate(
                    "",
                    xy=(row["end_pos"], pos),
                    arrowprops=dict(
                        facecolor="black", arrowstyle="->", connectionstyle="arc3", lw=2
                    ),
                )
            pos = pos + (np.max(z) / len(roi_df[roi_df["contig"] == i]))

        prophage_out_name = "hafeZ_prophages_for_" + str(i) + ".png"
        prophage_out_file: Path = Path(output) / prophage_out_name
        plt.savefig(prophage_out_file, format="png")
        plt.clf()


# define result
Result = collections.namedtuple("Result", ["protein", "hit", "bitscore", "evalue"])


def output_all_phrogs(
    roi_orf_aa: Dict[str, List[SeqRecord]],
    best_phrog_results_dict: Dict[str, Result],
    output: Path,
    database: Path,
    evalue: float,
) -> None:
    """
    Output the results of phrog annotation for ROIs and ORFs.

    Args:
        roi_orf_aa (Dict[str, List[SeqRecord]]): A dictionary containing ROIs and associated ORF amino acid sequences.
        best_phrog_results_dict (Dict[str, Result]): A dictionary containing the best phrog results for each ORF.
        output (Path): The output directory.
        database (str): Path to the phrog database.
        evalue (float): E-value threshold for annotation.

    Returns:
        None
    """

    # Initialize an empty DataFrame
    df = pd.DataFrame(columns=["roi", "orf", "evalue"])

    # Iterate through the best_phrog_results_dict and add rows to the DataFrame
    for roi, orf_records in roi_orf_aa.items():
        for orf_record in orf_records:
            new_row = {
                "roi": roi,
                "orf": orf_record.id,
            }
            df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)

    df["phrog"] = ""
    for index, row in df.iterrows():
        orf = row["orf"]
        if orf in best_phrog_results_dict.keys():  # where annotated
            orf_info = best_phrog_results_dict.get(orf, {})
            if orf_info.evalue < evalue:  # below eval threshold
                phrog_value = orf_info.hit
                eval = "{:.3e}".format(orf_info.evalue)
                df.at[index, "phrog"] = phrog_value
                df.at[index, "evalue"] = eval
            else:
                df.at[index, "phrog"] = "No_Hit"
                df.at[index, "evalue"] = np.nan
        else:
            df.at[index, "phrog"] = "No_Hit"
            df.at[index, "evalue"] = np.nan

    # merge in annots

    annot_file: Path = Path(database) / "phrog_annot_v4.tsv"
    phrogs_annot_df = pd.read_csv(
        annot_file,
        sep="\t",
        usecols=[0, 2, 3],
        names=["phrog", "annot", "category"],
    )

    phrogs_annot_df["phrog"] = "phrog_" + phrogs_annot_df["phrog"].astype(str)

    df = df.merge(phrogs_annot_df, on="phrog", how="left")

    # Replace NaN values in the 'annot' column with 'hypothetical protein'
    df["annot"].fillna("hypothetical protein", inplace=True)

    # Replace NaN values in the 'category' column with 'unknown function'
    df["category"].fillna("unknown function", inplace=True)

    df = df.reset_index(drop=True)
    df.to_csv(output / "hafeZ_pyhmmer_hits.tsv", sep="\t", index=False)


def output_roi_orfs(
    roi_orf_dna: Dict[str, SeqRecord],
    roi_df: pd.DataFrame,
    output: Path,
    roi_orf_aa: Dict[str, SeqRecord],
) -> None:
    """
    Output ROI-associated ORFs in both DNA and amino acid sequences.

    Args:
        roi_orf_dna (Dict[str, SeqRecord]): A dictionary containing ROI-associated ORF DNA sequences.
        roi_df (pd.DataFrame): DataFrame containing ROI information.
        output (Path): The output directory.
        roi_orf_aa (Dict[str, SeqRecord]): A dictionary containing ROI-associated ORF amino acid sequences.

    Returns:
        None
    """

    for i in roi_df["roi"].unique():
        orf_list = []
        for j in roi_orf_dna[i]:
            name = roi_df["roi_new"][roi_df["roi"] == i].iloc[0]
            j.id = name + "~" + j.id.split("~")[-1]
            orf_list.append(j)

        dna_file: Path = Path(output) / ("hafeZ_orfs_dna_" + str(name) + ".fasta")
        SeqIO.write(orf_list, dna_file, "fasta")
    for i in roi_df["roi"].unique():
        orf_list = []
        for j in roi_orf_aa[i]:
            name = roi_df["roi_new"][roi_df["roi"] == i].iloc[0]
            j.id = name + "~" + j.id.split("~")[-1]
            orf_list.append(j)

        aa_file: Path = Path(output) / ("hafeZ_orfs_aa_" + str(name) + ".fasta")
        SeqIO.write(orf_list, aa_file, "fasta")


def output_contig_Z(depths, output, median, mad):
    matplotlib.use("Agg")
    
    for i in depths:
        z = [((0.6745 * (x - median)) / mad) for x in depths[i]]
        fig = plt.gcf()
        fig.set_size_inches(18.5, 10.5)
        plt.plot(z)

        plot_file: Path = Path(output) / ("zscores_for_contig" + str(i) + ".png")

        plt.savefig(plot_file, format="png")
        plt.clf()


def output_roi_table(roi_df: pd.DataFrame, output: Path, depths: dict) -> None:
    """
    Generate and save a summary table for regions of interest (ROIs).

    Args:
        roi_df (pd.DataFrame): A DataFrame containing ROI data.
        output (Path): The directory where the summary table will be saved.
        depths (dict): Dictionary of coverage depths for each contig.

    Returns:
        None
    """

    for index, row in roi_df.iterrows():
        i = row['contig']
        if row["circular"] == False:
            roi_df.loc[index, "roi_length"] = row["end_pos"] - row["start_pos"] + 1
            depth = np.mean(depths[i][(row["start_pos"]-1):(row["end_pos"] - 1)])
        elif (row["circular"] == True) & (row["roi"].split("_")[-1] == "c1"):
            if row["start_pos"] > row["end_pos"]:
                len_1 = row["contig_len"] - row["start_pos"]
                len_2 = row["end_pos"] - 0
                depth = np.mean(np.concatenate([depths[i][row["start_pos"]-1:(row["contig_len"]-1)], depths[i][0:row["end_pos"]-1]]))
                roi_df.loc[index, "roi_length"] = len_1 + len_2
            else:
                depth = np.mean(depths[i][(row["start_pos"]-1):(row["end_pos"] - 1)])
                roi_df.loc[index, "roi_length"] = row["end_pos"] - row["start_pos"] + 1
        elif (row["circular"] == True) & (row["roi"].split("_")[-1] == "c2"):
            if row["start_pos"] > row["end_pos"]:
                roi_df.loc[index, "roi_length"] = row["end_pos"] - row["start_pos"] + 1
                depth = np.mean(np.concatenate([depths[i][row["start_pos"]-1:(row["contig_len"]-1)], depths[i][0:row["end_pos"]-1]]))
            else:
                len_1 = row["contig_len"] - row["start_pos"]
                len_2 = row["end_pos"] - 0
                roi_df.loc[index, "roi_length"] = len_1 + len_2
                depth = np.mean(depths[i][(row["start_pos"]-1):(row["end_pos"] - 1)])
        
        # set depth
        roi_df.loc[index, "mean_depth"] = round(depth, 2) 

    roi_df = roi_df[
        [
            "roi_new",
            "start_pos",
            "start_count",
            "end_pos",
            "end_count",
            "roi_length",
            "orf_count",
            "frac_phrog",
            "circular",
            "attL_seq",
            "attR_seq",
            "mean_depth"
        ]
    ].copy()
    roi_df.columns = [
        "roi",
        "start_pos",
        "start_count",
        "end_pos",
        "end_count",
        "roi_length",
        "orf_count",
        "frac_phrog",
        "circular",
        "attL_seq",
        "attR_seq",
        "mean_depth"
    ]
    for index, row in roi_df.iterrows():
        roi_df.loc[index, "start_pos"] = row["start_pos"] + 1
        roi_df.loc[index, "frac_phrog"] = round(row["frac_phrog"], 2)
        roi_df.loc[index, "roi_length"] = int(row["roi_length"])

    num_rois = roi_df.shape[0]

    logger.info(
        f"hafeZ has found {num_rois} ROI(s) indicating likely induced prophages."
    )

    out_file: Path = Path(output) / "hafeZ_summary_all_rois.tsv"

    roi_df.to_csv(out_file, sep="\t", index=False)
