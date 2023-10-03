import time
from Bio import SeqIO
import numpy as np
import subprocess
import pandas as pd
from io import StringIO
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import re
from scipy import stats
import pysam
from collections import Counter
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import multiprocessing
from functools import partial
import itertools
import os
from loguru import logger
from typing import Dict, Union, Tuple, Any, List
from typing import Optional
from hafeZ.utils.util import split_df_for_multiproc
from hafeZ.utils.external_tools import ExternalTool
from hafeZ.utils.exit import exit_error_gracefully

class Haf:
    """Main Class for hafeZ RoI Processing"""

    def __init__(
        self,
        roi_df: pd.DataFrame() = pd.DataFrame(
            {"col1": [1, 2, 3], "col2": [4, 5, 6]}
        ),
        sam_df: pd.DataFrame() = pd.DataFrame(
            {"col1": [1, 2, 3], "col2": [4, 5, 6]}
        ),
        end_roi_df: pd.DataFrame() = pd.DataFrame(
            {"col1": [1, 2, 3], "col2": [4, 5, 6]}
        ),
        clips_df:  pd.DataFrame() = pd.DataFrame(
            {"col1": [1, 2, 3], "col2": [4, 5, 6]}
        ),
        circular_df:  pd.DataFrame() = pd.DataFrame(
            {"col1": [1, 2, 3], "col2": [4, 5, 6]}
        ),
        clip_roi_df:  pd.DataFrame() = pd.DataFrame(
            {"col1": [1, 2, 3], "col2": [4, 5, 6]}
        ),
        clip_end_df:  pd.DataFrame() = pd.DataFrame(
            {"col1": [1, 2, 3], "col2": [4, 5, 6]}
        ),
        zscores:  Dict[Union[str, int], List[float]] = {
            'key1': [1.0, 2.0, 3.0],
            2: [4.0, 5.0, 6.0],
            'key3': [7.0, 8.0, 9.0],
        },
        median = 0.0, 
        mad = 0.0, 
        all_zscores = False,
        output: Path = Path("hafez_out"),
        depths:  Dict[Union[str, int], List[float]] = {
            'key1': [1.0, 2.0, 3.0],
            2: [4.0, 5.0, 6.0],
            'key3': [7.0, 8.0, 9.0],
        },
        start_time = 0.0,
        sam_secondary_df: pd.DataFrame() = pd.DataFrame(
            {"col1": [1, 2, 3], "col2": [4, 5, 6]}
        ),
    ) -> None:
        """
        Parameters
        --------
        roi_df: pd.DataFrame
            dataframe containing ROI information
        sam_df: pd.DataFrame
            dataframe containing SAM file output
        end_roi_df: pd.DataFrame
            dataframe containing ROIs that span the whole contig
        clips_df: pd.DataFrame
            dataframe containing soft clip data
        circular_df: pd.DataFrame
            dataframe containing circular contigs
        clip_roi_df: pd.DataFrame
            dataframe containing processed clips for roi_df 
        clip_end_df: pd.DataFrame
            dataframe containing processed clips for roi_df 
        zscores: Dict[Union[str, int], List[float]]
            Dict containing zscores
        median: float
            median z score
        mad: float
            mad z score
        all_zscores: bool
            whether to save all z scores.
        output: Path
            output directory
        depths: Dict[Union[str, int], List[float]]
            Dict containing depths
        start_time: float 
            start time
        sam_secondary_df:  pd.DataFrame
            dataframe containing SAM file outputvfor secondary alignments
        """
        self.roi_df = roi_df
        self.sam_df = sam_df
        self.end_roi_df = end_roi_df
        self.clips_df = clips_df
        self.circular_df = circular_df
        self.clip_roi_df = clip_roi_df
        self.clip_end_df = clip_end_df
        self.zscores = zscores
        self.median = median
        self.mad = mad
        self.all_zscores = all_zscores
        self.output = output
        self.depths = depths
        self.start_time = start_time
        self.sam_secondary_df = sam_secondary_df






    def get_ROIs(self, cutoff: float, width: int) -> None:
        """
        Extract Regions of Interest (ROIs) from Z-score data.

        Parameters:
            zscores (Dict[Union[str, int], List[float]]): A dictionary containing Z-scores for each contig/region.
            cutoff (float): The Z-score threshold for defining ROIs.
            width (int): The minimum width (difference) between start and end positions of an ROI.

        Returns:
            None

        This method identifies Regions of Interest (ROIs) from Z-score data. It scans the Z-score data for each
        contig/region and identifies segments that rise above the specified 'cutoff' and meet the 'width' criteria.
        ROIs are recorded in a DataFrame containing information such as accession, start, end, length, contig_len,
        and contig.

        If no ROIs are found, an error message is logged.

        Example:
            zscores_data = {
                'contig1': [0.5, 0.8, 1.2, 1.0],
                'contig2': [0.3, 0.6, 1.5, 1.2]
            }
            self.get_ROIs(zscores_data, 0.7, 2)
        """
        x = {}
        roi_dict = {
            "accession": [],
            "start": [],
            "end": [],
            "length": [],
            "contig_len": [],
            "contig": [],
        }

        for i in self.zscores:
            x[i] = {"start": [], "end": []}
            above_threshold = np.diff(self.zscores[i] > cutoff, prepend=False)
            rise_above_thres = np.argwhere(above_threshold)[::2, 0]
            sink_below_thres = np.argwhere(above_threshold)[1::2, 0]
            contig_len = len(np.array(self.zscores[i]))
            x[i]["start"] = rise_above_thres
            x[i]["end"] = sink_below_thres

            for j in np.arange(0, len(x[i]["start"])):
                if x[i]["end"][j] - x[i]["start"][j] > width:
                    roi_dict["length"].append(x[i]["end"][j] - x[i]["start"][j])
                    roi_dict["accession"].append(
                        "{}~{}~{}~{}".format(i, x[i]["start"][j], x[i]["end"][j], j)
                    )
                    roi_dict["start"].append(x[i]["start"][j])
                    roi_dict["end"].append(x[i]["end"][j])
                    roi_dict["contig_len"].append(contig_len)
                    roi_dict["contig"].append(i.split("~")[0])

        self.roi_df = pd.DataFrame(roi_dict) if len(roi_dict) > 0 else None

        # exit 
        if self.roi_df  is  None:
            exit_error_gracefully(self.output, self.all_zscores, self.depths, self.median, self.mad, self.start_time)
            logger.error("No RoIs found at all.")
        

    def join_ROIs(self, join_window: int) -> None:
        """
        Merge close Regions of Interest (ROIs).

        join_window (int): window where ROIs will be merged (10000 default). 

        This method merges ROIs that are close to each other on the same contig. It scans the ROIs in the
        'self.roi_df' DataFrame and combines ROIs that are within 10,000 base pairs of each other along the
        same contig. The merged ROIs are then stored in 'self.roi_df' and returned as a DataFrame.

        Returns:
           None
        """

        counter_1 = 0
        joined_roi_df = {
            "accession": [],
            "start": [],
            "end": [],
            "length": [],
            "contig_len": [],
            "contig": [],
        }
        for index, row in self.roi_df.iterrows():
            counter_1 = counter_1 + 1
            past_df_len = 0
            sub_df = self.roi_df[
                (self.roi_df["contig"] == row["contig"])
                & (
                    (
                        (self.roi_df["start"].between(row["start"] - join_window, row["end"] + join_window))
                        | (self.roi_df["end"].between(row["start"] - join_window, row["end"] + join_window))
                    )
                    | (
                        (
                            (row["start"] >= self.roi_df["start"] - join_window)
                            & (row["start"] <= self.roi_df["end"] + join_window)
                        )
                        | (
                            (row["end"] >= self.roi_df["start"] - join_window)
                            & (row["end"] <= self.roi_df["end"] + join_window)
                        )
                    )
                )
            ]
            current_df_len = len(sub_df)
            sub_list = list(sub_df["start"]) + list(sub_df["end"])
            min_pos = min(sub_list)
            max_pos = max(sub_list)
            while current_df_len > past_df_len:
                past_df_len = current_df_len
                sub_df = self.roi_df[
                    (self.roi_df["contig"] == row["contig"])
                    & (
                        (
                            self.roi_df["start"].between(min_pos - join_window, max_pos + join_window)
                            | self.roi_df["end"].between(min_pos - join_window, max_pos + join_window)
                        )
                        | (
                            (
                                (min_pos >= self.roi_df["start"] - join_window)
                                & (min_pos <= self.roi_df["end"] + join_window)
                            )
                            | (
                                (max_pos >= self.roi_df["start"] - join_window)
                                & (max_pos <= self.roi_df["end"] + join_window)
                            )
                        )
                    )
                ]
                sub_list = list(sub_df["start"]) + list(sub_df["end"])
                min_pos = min(sub_list)
                max_pos = max(sub_list)
                current_df_len = len(sub_df)
            joined_roi_df["accession"].append(
                "~".join(sub_df["accession"].iloc[0].split("~")[:-1]) + "~" + str(counter_1)
            )
            joined_roi_df["start"].append(min_pos)
            joined_roi_df["end"].append(max_pos)
            joined_roi_df["length"].append(max_pos - min_pos)
            joined_roi_df["contig_len"].append(sub_df["contig_len"].iloc[0])
            joined_roi_df["contig"].append(sub_df["contig"].iloc[0])

        joined_roi_df = pd.DataFrame.from_dict(joined_roi_df)
        joined_roi_df = joined_roi_df.drop_duplicates(
            subset=["start", "end", "length", "contig_len", "contig"], keep="first"
        )
        self.roi_df = joined_roi_df


    def filter_ROIs(self, width: int) -> None:
        """
        Filter small Regions of Interest (ROIs).

        Parameters:
            width (int): The minimum width (difference) for ROIs to be retained.

        Returns:
            None

        This method filters ROIs in the 'self.roi_df' DataFrame based on their length (width). ROIs with a length
        less than 'width' are removed from the DataFrame. The filtered ROIs are then returned as a new DataFrame.

        Example:
            self.filter_ROIs(1000)
        """


        filtered_roi_df = self.roi_df[self.roi_df["length"] >= width]

        if len(filtered_roi_df) == 0: 
            exit_error_gracefully(self.output, self.all_zscores, self.depths, self.median, self.mad, self.start_time)
            logger.error("After discarding small ROIs, no ROIs found.")

        self.roi_df = filtered_roi_df



    def bamfile_to_pandas(self, output: Path, neighbourhood: int) -> None:
        """
        Parse and process a sorted BAM file into a Pandas DataFrame.

        Parameters:
            output (Path): The output folder where the SAM file is located.
            neighbourhood (int): number of bases before and after the ROI
        
        Returns:
            None

        This method parses and processes a sorted BAM file into a Pandas DataFrame. It filters reads based on the ROIs provided
        in 'self.roi_df', ensuring that only reads within the specified regions are included in the output DataFrame.
        
        """

        sorted_bam: Path = output / "temp_minimap_sorted.bam"

        bamfile = pysam.AlignmentFile(sorted_bam, "rb")

        df_dict: Dict[str, list] = {
            "qname": [],
            "rname": [],
            "pos": [],
            "cigar": [],
            "tlen": [],
            "seq": [],
            "flag": [],
            "pnext": [],
            "qqual": [],
        }

        for index, row in self.roi_df.iterrows():
            contig = row["accession"].split("~")[0]
            start = row["start"] - neighbourhood
            if start < 0:
                start = 0
            end = row["end"] + neighbourhood
            if end > row["contig_len"]:
                end = row["contig_len"]
            for reads in bamfile.fetch(contig, start, end):
                df_dict["qname"].append(reads.query_name)
                df_dict["rname"].append(contig)
                df_dict["pos"].append(reads.reference_start)
                df_dict["cigar"].append(reads.cigarstring)
                df_dict["tlen"].append(reads.template_length)
                df_dict["seq"].append(reads.query_sequence)
                df_dict["flag"].append(reads.flag)
                df_dict["pnext"].append(reads.next_reference_start)
                df_dict["qqual"].append(reads.query_qualities.tolist())

        df = pd.DataFrame.from_dict(df_dict)
        df = df.drop_duplicates(
            subset=["qname", "rname", "pos", "cigar", "tlen", "seq", "flag", "pnext"],
            keep="first",
        )
        df["tlen"] = df["tlen"].astype(int)
        df["pos"] = df["pos"].astype(int)

        if df.empty:
            exit_error_gracefully(self.output, self.all_zscores, self.depths, self.median, self.mad, self.start_time)
            logger.error("No depth found at all positions.")

        self.sam_df = df



    def find_contig_end_rois(self, neighbourhood: int) -> Tuple[Union[pd.DataFrame, str], Union[pd.DataFrame, str]]:
        """
        Find ROIs at the ends of contigs and categorize them.

        Returns:
            Tuple[Union[pd.DataFrame, str], Union[pd.DataFrame, str]]: A tuple containing two elements:
                - DataFrame containing ROIs at the ends of contigs and their categorization.
                - Updated DataFrame of ROIs after removing those at the ends of contigs.
            neighbourhood (int): number of bases before and after the ROI

        This method identifies ROIs at the ends of contigs in 'self.roi_df', categorizes them, and separates them into
        a separate DataFrame 'end_df'. The remaining ROIs are updated in 'roi_df'. The method returns both DataFrames
        or "exit" if there are no ROIs.

        Example:
            end_rois, remaining_rois = self.find_contig_end_rois(10000)
        """
        for index, row in self.roi_df.iterrows():
            self.roi_df.loc[index, "center"] = (row["start"] + row["end"]) / 2

        end_df = self.roi_df[
            (
                (self.roi_df["start"].between(0, neighbourhood))
                | (
                    self.roi_df["end"].between(
                        self.roi_df["contig_len"] - neighbourhood, self.roi_df["contig_len"]
                    )
                )
            )
        ].copy()

        self.roi_df = self.roi_df[
            ~(
                (self.roi_df["start"].between(0, neighbourhood))
                | (
                    self.roi_df["end"].between(
                        self.roi_df["contig_len"] - neighbourhood, self.roi_df["contig_len"]
                    )
                )
            )
        ].copy()

        for index, row in end_df.iterrows():
            if (0 <= row["start"] <= neighbourhood) and (
                row["contig_len"] - neighbourhood <= row["end"] <= row["contig_len"]
            ):
                end_df.loc[index, "end_info"] = "whole"
            elif (0 <= row["start"] <= neighbourhood) and (
                not row["contig_len"] - neighbourhood <= row["end"] <= row["contig_len"]
            ):
                end_df.loc[index, "end_info"] = "left"
            elif (not 0 <= row["start"] <= neighbourhood) and (
                row["contig_len"] - neighbourhood <= row["end"] <= row["contig_len"]
            ):
                end_df.loc[index, "end_info"] = "right"

        self.roi_df["end_info"] = np.nan

        if len(self.roi_df) < 1:
            exit_error_gracefully(self.output, self.all_zscores, self.depths, self.median, self.mad, self.start_time)
            logger.error("Error with finding ROIs that constitute a whole contig.")

        if len(end_df) < 1:
            logger.info("Zero ROIs found that span whole contigs.")

        # end_df saved to the class
        self.end_roi_df = end_df


    def find_far_reads(self, width: int, neighbourhood: int) -> Union[pd.DataFrame, str]:
        """
        Find reads that are far from ROIs and update the ROI DataFrame.

        Parameters:
            width (int): The width threshold to consider reads as far from ROIs.
            neighbourhood (int): number of bases before and after the ROI

        Returns:
            Union[pd.DataFrame, str]: If ROIs are found, returns the updated DataFrame of ROIs.
            If no ROIs are found, logger errors.

        This method identifies reads that are far from ROIs based on the specified width threshold. It updates the
        'self.roi_df' DataFrame to include information about far reads, such as their positions and counts.

        Example:
            updated_roi_df = self.find_far_reads(500, 15001)
        """

        df = self.sam_df[
            ((self.sam_df["tlen"] > width) | (self.sam_df["tlen"] < (width * -1)))
            & (self.sam_df.duplicated(subset=["qname"], keep=False))
        ].copy()

        for index, row in self.roi_df.iterrows():
            center = row["center"]
            sub_df = df[
                (df["rname"] == row["contig"])
                & (
                    (
                        (df["pos"].between(row["start"] - neighbourhood, center))
                        & (df["pnext"].between(center, row["end"] + neighbourhood))
                    )
                    | (
                        (df["pnext"].between(row["start"] - neighbourhood, center))
                        & (df["pos"].between(center, row["end"] + neighbourhood))
                    )
                )
            ].copy()

            sub_df = sub_df[sub_df.duplicated(subset=["qname"], keep=False)]

            if (len(sub_df) / 2) < 1:
                self.roi_df = self.roi_df[self.roi_df["accession"] != row["accession"]].copy()
            else:
                self.roi_df.loc[index, "far_start"] = np.median(
                    sub_df["pos"][sub_df["pos"] < center].copy()
                )
                self.roi_df.loc[index, "far_end"] = np.median(
                    sub_df["pos"][sub_df["pos"] > center].copy()
                )
                self.roi_df.loc[index, "far_count_pairs"] = len(sub_df.copy()) / 2

        if len(self.roi_df) < 1:
            logger.info("No ROIs with evidence of pairs of reads that are distant from each other.")


    def find_contig_end_rois_again(self, neighbourhood: int) -> None:
        """
        Check for ROIs near the ends of contigs again and update the ROI DataFrame.

        parameters:
            neighbourhood (int): number of bases before and after the ROI

        Returns:
            None

        This method checks for ROIs near the ends of contigs in the self.roi_df DataFrame and updates
        the 'end_info' column accordingly. It then combines the updated ROIs with the original DataFrame.

        """

        for index, row in self.roi_df.iterrows():
            self.roi_df.loc[index, "center"] = (row["start_pos"] + row["end_pos"]) / 2

        end_df = self.roi_df[
            (
                (self.roi_df["start_pos"].between(0, neighbourhood))
                | (
                    self.roi_df["end_pos"].between(
                        self.roi_df["contig_len"] - neighbourhood, self.roi_df["contig_len"]
                    )
                )
            )
        ].copy()

        self.roi_df = self.roi_df[
            ~(
                (self.roi_df["start_pos"].between(0, neighbourhood))
                | (
                    self.roi_df["end_pos"].between(
                        self.roi_df["contig_len"] - neighbourhood, self.roi_df["contig_len"]
                    )
                )
            )
        ].copy()

        for index, row in end_df.iterrows():
            if (0 <= row["start_pos"] <= neighbourhood) and (
                row["contig_len"] - neighbourhood <= row["end_pos"] <= row["contig_len"]
            ):
                end_df.loc[index, "end_info"] = "whole"
            elif (0 <= row["start_pos"] <= neighbourhood) and (
                not row["contig_len"] - neighbourhood <= row["end_pos"] <= row["contig_len"]
            ):
                end_df.loc[index, "end_info"] = "left"
            elif (not 0 <= row["start_pos"] <= neighbourhood) and (
                row["contig_len"] - neighbourhood <= row["end_pos"] <= row["contig_len"]
            ):
                end_df.loc[index, "end_info"] = "right"

        self.roi_df["end_info"] = np.nan
        self.roi_df = pd.concat([self.roi_df, end_df])
        self.roi_df = self.roi_df.reset_index(drop=True)

    def find_soft_clippings(
        self, width: int, threads: int, neighbourhood: int
    ) -> Tuple[Optional[pd.DataFrame], pd.DataFrame]:
        """
        Find soft-clipped reads within ROIs and create a DataFrame with the results.

        Parameters:
            width (int): The minimum width threshold for soft clippings to be considered.
            threads (int): Number of threads for parallel processing.
            neighbourhood (int): Number of bases in the neighbourhood

        Returns:
            Tuple[Optional[pd.DataFrame], pd.DataFrame]: A tuple containing the DataFrame with soft clippings
            meeting the width threshold and the DataFrame of soft clippings (clips). Returns None for the first
            element if no soft clippings are found.

        This method identifies soft-clipped reads within the provided ROIs and creates a DataFrame
        with information about their positions and counts. It also returns the DataFrame of soft clippings (clips).
        The function uses parallel processing with the specified number of threads for efficient computation.

        Example:
            soft_clippings_df, clips_df = self.find_soft_clippings(20, 4)
        """

        clips = self.sam_df[
            (self.sam_df["cigar"].str.count("M") == 1)
            & (self.sam_df["cigar"].str.count("S") == 1)
            & (self.sam_df["cigar"].str.count("I") == 0)
            & (self.sam_df["cigar"].str.count("X") == 0)
            & (self.sam_df["cigar"].str.count("D") == 0)
        ]

        chunks, pool = split_df_for_multiproc(self.roi_df, threads)

        df = pd.concat(
            pool.map(
                partial(multiprocess_find_soft_clippings, clips=clips, width=width, neighbourhood=neighbourhood),
                chunks,
            )
        )

        df = df.reset_index(drop=True)
        
        self.roi_df = df
        self.clips_df = clips



    def calc_roi_Z_medians(
        self,
        cov,
        median,
        mad,
        median_z_cutoff,
        threads
    ) -> pd.DataFrame:
        """
        Calculate the median Z-scores for each region of interest (ROI).

        Args:
            cov (dict): A dictionary containing coverage information for each contig.
            median (float): The median depth value for all contigs.
            mad (float): The median absolute deviation (MAD) for all contigs.
            median_z_cutoff (float): The Z-score cutoff value for considering a region.
            threads (int): The number of threads to use for parallel processing.

        Returns:
            pd.DataFrame: DataFrame with calculated median Z-scores for each ROI.

        This method calculates the median Z-scores for each region of interest (ROI) in the
        self.roi_df DataFrame based on the provided depth information, median, MAD, and Z-score cutoff.
        It returns a DataFrame with the calculated Z-scores.

        """

        self.roi_df = self.roi_df.reset_index(drop=True)
        self.roi_df["start_pos"] = self.roi_df["start_pos"].astype("int")
        self.roi_df["end_pos"] = self.roi_df["end_pos"].astype("int")

        chunks, pool = split_df_for_multiproc(self.roi_df, threads)
        self.roi_df = pd.concat(
            pool.map(
                partial(
                    multiprocessing_calc_roi_z_medians,
                    depths=cov,
                    median=median,
                    mad=mad,
                    median_z_cutoff=median_z_cutoff,
                ),
                chunks,
            )
        )


    def filter_Z_medians(self) -> None:
        """
        Filter regions of interest (ROIs) by Z-score.

        Returns:
            None

        This method filters ROIs in self.roi_df based on the provided Z-score cutoff.
        It returns a DataFrame with the filtered ROIs.

        """

        both_ends_df = self.roi_df[
            (self.roi_df["end_info"] == "whole") & (self.roi_df["longest_below_z"] < 7500)
        ]
        left_df = self.roi_df[
            (self.roi_df["end_info"] == "left")
            & (self.roi_df["longest_below_z"] < 7500)
            & (
                (self.roi_df["med_z"] >= self.roi_df["med_z"] / 2)
                & (self.roi_df["med_z_r"] < self.roi_df["med_z"] / 2)
                & (self.roi_df["med_z_l"].isna())
                | (self.roi_df["med_z"] >= self.roi_df["med_z"] / 2)
                & (self.roi_df["med_z_r"] < self.roi_df["med_z"] / 2)
                & (self.roi_df["med_z_l"] < self.roi_df["med_z"] / 2)
            )
        ]
        right_df = self.roi_df[
            (self.roi_df["end_info"] == "right")
            & (self.roi_df["longest_below_z"] < 7500)
            & (
                (self.roi_df["med_z"] >= self.roi_df["med_z"] / 2)
                & (self.roi_df["med_z_l"] < self.roi_df["med_z"] / 2)
                & (self.roi_df["med_z_r"].isna())
                | (self.roi_df["med_z"] >= self.roi_df["med_z"] / 2)
                & (self.roi_df["med_z_r"] < self.roi_df["med_z"] / 2)
                & (self.roi_df["med_z_l"] < self.roi_df["med_z"] / 2)
            )
        ]
        contained_df = self.roi_df[
            (self.roi_df["med_z"] >= self.roi_df["med_z"] / 2)
            & (self.roi_df["longest_below_z"] < 7500)
            & (self.roi_df["med_z_l"] < self.roi_df["med_z"] / 2)
            & (self.roi_df["med_z_r"] < self.roi_df["med_z"] / 2)
            & (self.roi_df["end_info"].isna())
        ]

        self.roi_df = pd.concat([both_ends_df, left_df, right_df, contained_df])
        self.roi_df = self.roi_df.reset_index(drop=True)

        if len(self.roi_df) < 1:
            exit_error_gracefully(self.output, self.all_zscores, self.depths, self.median, self.mad, self.start_time)
            logger.error("No ROIs found after Z-score filtering.")


    def keep_only_x_best(self, keep_threshold):
        """
        Keep only the top "best" ROI possibilities.

        Args:
            keep_threshold (int): The maximum number of ROIs to keep.

        Returns:
            pd.DataFrame: DataFrame with the top "best" ROIs.

        This method filters and keeps only the top "best" ROI possibilities based on
        the specified keep_threshold. It returns a DataFrame with the selected ROIs.

        """

        df_list = []

        for i in self.roi_df["roi"].unique():
            sub_df = self.roi_df[self.roi_df["roi"] == i]
            sub_df = sub_df.reset_index(drop=True)

            for index, row in sub_df.iterrows():
                sub_df.loc[index, "total_clips"] = int(row["start_count"]) + int(
                    row["end_count"]
                )

            sub_df = sub_df.sort_values(by=["total_clips"], ascending=False).iloc[
                0:keep_threshold, :
            ]

            df_list.append(sub_df)

        self.roi_df = pd.concat(df_list)
        self.roi_df = self.roi_df.reset_index(drop=True)



    def check_for_circular_roi(self) -> None:
        """
        Check for circular regions of interest (ROIs) in the context of plasmid phages.

        Args:
            sam_df (pd.DataFrame): DataFrame containing SAM file data.

        Returns:
            If no circular ROIs are found, the circular ROI DataFrame will be None
            If no non-circular ROIs are left, the non-circular ROI DataFrame will be None

        This method checks for circular regions of interest (ROIs) in the context of plasmid phages.
        It returns a tuple of two DataFrames, one containing circular ROIs and the other containing non-circular ROIs.
        If no circular ROIs are found, the circular ROI DataFrame will be None
        If no non-circular ROIs are left, the non-circular ROI DataFrame will be None

        Example:
            circular_df, non_circular_df = self.check_for_circular_roi()
        """
        logger.info(" Checking for plasmid phages.")
        overlaps_dict = {}
        self.roi_df["circular"] = False

        for index, row in self.roi_df[~pd.isna(self.roi_df["end_info"])].iterrows():
            start = row["start_pos"]
            end = row["end_pos"]
            contig_len = row["contig_len"]
            contig = row["roi"].split("~")[0]
            read_list = []

            if start < 1000:
                read_list = list(
                    self.sam_df["qname"][
                        (self.sam_df["rname"] == contig) & (self.sam_df["pos"].between(0, 1000))
                    ]
                )
            elif end > contig_len - 1000:
                read_list = list(
                    self.sam_df["qname"][
                        (self.sam_df["rname"] == contig)
                        & (self.sam_df["pos"].between(contig_len - 1000, contig_len))
                    ]
                )

            sub_df = self.sam_df[
                (self.sam_df["qname"].isin(read_list)) & (self.sam_df["rname"] == contig)
            ]

            sub_df = sub_df[
                (
                    (sub_df["pos"].between(0, 1000))
                    & (sub_df["pnext"].between(contig_len - 1000, contig_len))
                )
                | (
                    (sub_df["pnext"].between(0, 1000))
                    & (sub_df["pos"].between(contig_len - 1000, contig_len))
                )
            ]

            sub_df = sub_df[sub_df.duplicated(subset=["qname", "rname"], keep=False)]

            min_pos = sub_df["pos"][sub_df["pos"] <= 1000].max()
            max_pos = sub_df["pos"][sub_df["pos"] >= contig_len - 1000].min()

            overlaps = self.roi_df["roi"][
                (self.roi_df["roi"].str.split("~").str[0] == contig)
                & (
                    (self.roi_df["start_pos"].between(0, min_pos))
                    | (self.roi_df["end_pos"].between(max_pos, contig_len))
                )
            ].unique()

            if len(overlaps) >= 1:
                overlaps_dict[row["roi"]] = overlaps

        overlaps = []
        for i in overlaps_dict:
            overlaps.append([k for (k, v) in overlaps_dict.items() if i in v])

        unique_overlaps = [list(x) for x in set(tuple(x) for x in overlaps)]
        counter = 1000

        for overlaps in unique_overlaps:
            starts = []
            ends = []

            for i in overlaps:
                starts.append(self.roi_df["start_pos"][self.roi_df["roi"] == i].iloc[0])
                ends.append(self.roi_df["end_pos"][self.roi_df["roi"] == i].iloc[0])
                contig_len = self.roi_df["contig_len"][self.roi_df["roi"] == i].iloc[0]
                self.roi_df = self.roi_df[self.roi_df["roi"] != i]

            contig = i.split("~")[0]
            start = np.max(starts)
            end = np.min(ends)
            roi_1 = contig + "~roi_" + str(counter) + "_c1"
            roi_2 = contig + "~roi_" + str(counter) + "_c2"
            df2 = {
                "roi": [roi_1, roi_2],
                "start_pos": [start, end],
                "start_count": [np.nan, np.nan],
                "end_pos": [end, start],
                "end_count": [np.nan, np.nan],
                "contig_len": [contig_len, contig_len],
                "contig": [contig, contig],
                "circular": [True, True],
            }
            df2 = pd.DataFrame.from_dict(df2)
            self.roi_df = pd.concat([self.roi_df, df2])
            counter = counter - 1

        self.circular_df = self.roi_df[self.roi_df["circular"] == True]
        self.roi_df = self.roi_df[self.roi_df["circular"] == False]

        if len(self.roi_df) < 1:
            exit_error_gracefully(self.output, self.all_zscores, self.depths, self.median, self.mad, self.start_time)
            logger.error("Error in finding circular ROIs (possible phage plasmids).")
        if len(self.circular_df) < 1:
            self.circular_df = None
            logger.info("No circular ROIs (possible phage plasmids) found.")
        else:
            for index, row in self.circular_df.iterrows():
                self.circular_df.loc[index, "contig_split"] = np.nan


    def collecting_clipped_reads(
        self,
        output: Path,
        threads: int
    ) -> None:
        """
        Collect soft clipped reads and write them to a FASTQ file.

        Args:
            clips (pd.DataFrame): DataFrame containing clipped reads data.
            output (Path): Path to the output folder where the FASTQ file will be saved.
            threads (int): Number of threads for parallel processing.

        Returns:
            None

        This method collects soft clipped reads and writes them to a FASTQ file.
        The collected reads are also removed from the ROI DataFrame.

        Example:
            updated_roi_df = self.collecting_clipped_reads(clips, output, threads)
        """
        self.roi_df = self.roi_df.reset_index(drop=True)
        chunks, pool = split_df_for_multiproc(self.roi_df, threads)
        values = pool.map(partial(collect_clips, clips=self.clips_df), chunks)
        read_list = list(itertools.chain.from_iterable([item[0] for item in values]))
        self.roi_df = pd.concat([item[1] for item in values])
        clipped_reads: Path = output / "temp_soft_clipped.fastq"
        SeqIO.write(read_list, clipped_reads, "fastq")
        


    def map_clipped_reads(self, output, threads, genome, memory_limit, logdir) -> None:
        """
        Map soft clipped reads to a reference genome using minimap2.

        Args:
            output (Path): Path to the output directory.
            threads (int): Number of threads for parallel processing.
            genome (str): Path to the reference genome file.
            memory_limit (str): Memory limit for samtools sort.
            logdir (Path): Path to the directory for log files.

        This method maps soft clipped reads to a reference genome using minimap2.
        It then sorts and indexes the resulting SAM file using samtools.

        """

        clipped_reads: Path = output / "temp_soft_clipped.fastq"
        out_sam: Path = output / "temp_minimap_clipped.sam"
        out_sorted_bam: Path = output / "temp_minimap_clipped_sorted.bam"

        minimap2 = ExternalTool(
        tool="minimap2",
        input=f"",
        output=f"",
        params=f'-ax sr {genome} {clipped_reads} -t {str(threads)} -o {out_sam} ',
        logdir=logdir,
    )
        
        ExternalTool.run_tool(minimap2)


        samtools_sort = ExternalTool(
        tool="samtools",
        input=f"",
        output=f"",
        params=f'sort  -@ {str(threads)} -m {memory_limit} -o {out_sorted_bam} {out_sam} ',
        logdir=logdir,
    )
        
        ExternalTool.run_tool(samtools_sort)

        samtools_index = ExternalTool(
        tool="samtools",
        input=f"",
        output=f"",
        params=f'index  -@ {str(threads)} -b {out_sorted_bam} ',
        logdir=logdir,
    )
        
        ExternalTool.run_tool(samtools_index)



    def seperate_rois_near_ends(self) -> None:
        """
        Separate ROIs near the ends of contigs.

        Returns:
            None

        This method separates ROIs that are near the ends of contigs into a separate DataFrame.
        It checks if ROIs span the entire contig or are close to the contig ends.
        If no qualifying ROIs are found, the corresponding DataFrames are set to None.

        """

        self.end_roi_df = self.roi_df[
            (
                (self.roi_df["start_pos"].between(0, 2000))
                | (
                    self.roi_df["end_pos"].between(
                        self.roi_df["contig_len"] - 2000, self.roi_df["contig_len"] + 1
                    )
                )
            )
        ].copy()  # pull out ROIs into a separate DataFrame that would have failed later tests unfairly due to being close to the end of the contig
        self.roi_df = self.roi_df[
            ~(
                (
                    (self.roi_df["start_pos"].between(0, 2000))
                    | (
                        self.roi_df["end_pos"].between(
                            self.roi_df["contig_len"] - 2000, self.roi_df["contig_len"] + 1
                        )
                    )
                )
            )
        ].copy()
        for index, row in self.end_roi_df.iterrows():
            if (0 <= row["start_pos"] <= 2000) and (
                row["contig_len"] - 2000 <= row["end_pos"] <= row["contig_len"]
            ):
                self.end_roi_df.loc[index, "end_info"] = "whole"
            elif (0 <= row["start_pos"] <= 2000) and (
                not row["contig_len"] - 2000 <= row["end_pos"] <= row["contig_len"]
            ):
                self.end_roi_df.loc[index, "end_info"] = "left"
            elif (not 0 <= row["start_pos"] <= 2000) and (
                row["contig_len"] - 2000 <= row["end_pos"] <= row["contig_len"]
            ):
                self.end_roi_df.loc[index, "end_info"] = "right"
        self.roi_df["contig_split"] = np.nan
        self.roi_df["end_info"] = np.nan
        if len(self.roi_df) < 1:
            self.roi_df = None
            exit_error_gracefully(self.output, self.all_zscores, self.depths, self.median, self.mad, self.start_time)
            logger.error("Error with finding ROIs near the end of contigs.No ROIs found.")
        if len(self.end_roi_df) < 1:
            self.end_roi_df = None
            


    def process_clip_sam(self, output: str) -> None:
        """
        Process the clipped SAM file to extract relevant information for ROIs.

        Args:
            output (Path): The folder containing the clipped SAM file.

        Returns:
            None
        """
        out_sorted_bam: Path = output / "temp_minimap_clipped_sorted.bam"
        samfile = pysam.AlignmentFile(
            out_sorted_bam, "rb"
        )
        df_dict = {
            "qname": [],
            "rname": [],
            "pos": [],
            "cigar": [],
            "tlen": [],
            "seq": [],
            "flag": [],
            "pnext": [],
        }
        for index, row in self.roi_df.iterrows():
            contig = row["contig"]
            start = row["start_pos"] - 200
            if start < 0:
                start = 0
            end = row["end_pos"] + 200
            if end > row["contig_len"]:
                end = row["contig_len"]
            for reads in samfile.fetch(contig, start, end):
                df_dict["qname"].append(reads.query_name)
                df_dict["rname"].append(contig)
                df_dict["pos"].append(reads.reference_start)
                df_dict["cigar"].append(reads.cigarstring)
                df_dict["tlen"].append(reads.template_length)
                df_dict["seq"].append(reads.query_sequence)
                df_dict["flag"].append(reads.flag)
                df_dict["pnext"].append(reads.next_reference_start)
        df = pd.DataFrame.from_dict(df_dict)
        df = df.drop_duplicates(
            subset=["qname", "rname", "pos", "cigar", "tlen", "seq", "flag", "pnext"],
            keep="first",
        )
        df["tlen"] = df["tlen"].astype(int)
        df["pos"] = df["pos"].astype(int)

        self.clip_roi_df = df


    def get_clip_pos(self, output: Path) -> None:
        """
        Get positions for clipped ends (contained ROIs) and calculate clip counts.

        Args:
            output (Path): The folder containing the output files.

        Returns:
            None
        """

        out_sorted_bam: Path = output / "temp_minimap_clipped_sorted.bam"

        samfile = pysam.AlignmentFile(
            out_sorted_bam, "rb"
        )
        names = pysam.IndexedReads(samfile)
        names.build()
        self.roi_df["start-end_clip_count"] = np.nan
        self.roi_df["end-start_clip_count"] = np.nan
        self.roi_df["total_clip_count"] = np.nan
        clip_df2 = self.clip_roi_df[
            (self.clip_roi_df["cigar"].str.count("M") == 1)
            & (self.clip_roi_df["cigar"].str.count("S") == 0)
            & (self.clip_roi_df["cigar"].str.count("I") == 0)
            & (self.clip_roi_df["cigar"].str.count("X") == 0)
            & (self.clip_roi_df["cigar"].str.count("D") == 0)
        ]
        for index, row in self.roi_df.iterrows():
            start_name = row["roi"] + "~" + str(row["start_pos"])
            end_name = row["roi"] + "~" + str(row["end_pos"])
            start_mapped_df = clip_df2[
                (clip_df2["qname"].str.split("__").str[0] == start_name)
                & (clip_df2["pos"].between(row["end_pos"] - 200, row["end_pos"] + 200))
            ]
            end_mapped_df = clip_df2[
                (clip_df2["qname"].str.split("__").str[0] == end_name)
                & (clip_df2["pos"].between(row["start_pos"] - 200, row["start_pos"] + 200))
            ]
            if (len(start_mapped_df) >= 10) and (len(end_mapped_df) >= 10):
                self.roi_df.loc[index, "start-end_clip_count"] = len(end_mapped_df)
                self.roi_df.loc[index, "end-start_clip_count"] = len(start_mapped_df)
                self.roi_df.loc[index, "total_clip_count"] = (
                    len(end_mapped_df)
                    + len(start_mapped_df)
                    + row["start_count"]
                    + row["end_count"]
                )  # added these last 2 to test something
            else:
                if (len(start_mapped_df) < 10) and (len(end_mapped_df) >= 10):
                    self.roi_df.loc[index, "start-end_clip_count"] = len(end_mapped_df)
                    self.roi_df.loc[index, "end-start_clip_count"] = check_read_lengths(
                        row, index, self.clip_roi_df, start_name, names, "start"
                    )
                    self.roi_df.loc[index, "total_clip_count"] = (
                        10 + len(start_mapped_df) + row["start_count"] + row["end_count"]
                    )  # added these last 2 to test something
                elif (len(end_mapped_df) < 10) and (len(start_mapped_df) >= 10):
                    self.roi_df.loc[index, "end-start_clip_count"] = len(start_mapped_df)
                    self.roi_df.loc[index, "start-end_clip_count"] = check_read_lengths(
                        row, index, self.clip_roi_df, end_name, names, "end"
                    )
                    self.roi_df.loc[index, "total_clip_count"] = (
                        len(end_mapped_df) + 10 + row["start_count"] + row["end_count"]
                    )  # added these last 2 to test something
        self.roi_df = self.roi_df[
            (self.roi_df["start-end_clip_count"].notna())
            & (self.roi_df["end-start_clip_count"].notna())
        ]
        if len(self.roi_df) < 1:
            exit_error_gracefully(self.output, self.all_zscores, self.depths, self.median, self.mad, self.start_time)
            logger.error("Error with finding ROIs with clipped reads. No ROIs found.")

        self.roi_df = self.roi_df.sort_values(by=["roi", "total_clip_count"], ascending=False)
        df_list = []
        for i in self.roi_df["roi"].unique():
            df_list.append(self.roi_df[self.roi_df["roi"] == i].iloc[[0]])
        self.roi_df = pd.concat(df_list)




    def process_clip_sam_end_rois(self, output: str):
        """
        Get positions for clipped ends (contained ROIs) and calculate clip counts for ends.

        Args:
            clip_df (pd.DataFrame): DataFrame containing clipped reads.
            output (Path): The folder containing the output files.

        Returns:
            None
        """
        out_sorted_bam: Path = output / "temp_minimap_clipped_sorted.bam"
        samfile = pysam.AlignmentFile(
            out_sorted_bam, "rb"
        )
        names = pysam.IndexedReads(samfile)
        names.build()
        df_dict = {
            "qname": [],
            "rname": [],
            "pos": [],
            "cigar": [],
            "tlen": [],
            "seq": [],
            "flag": [],
            "pnext": [],
            "roi": [],
        }
        for index, row in self.roi_df.iterrows():
            start_list = np.core.defchararray.add(
                "{}~{}__".format(row["roi"], str(row["start_pos"])),
                np.arange(1, row["start_count"] + 1, 1).astype(int).astype(str),
            )
            for read in start_list:
                try:
                    start_reads = names.find(read)
                    for x in start_reads:
                        df_dict["qname"].append(x.query_name)
                        df_dict["rname"].append(x.reference_name)
                        df_dict["pos"].append(x.reference_start)
                        df_dict["cigar"].append(x.cigarstring)
                        df_dict["tlen"].append(x.template_length)
                        df_dict["seq"].append(x.query_sequence)
                        df_dict["flag"].append(x.flag)
                        df_dict["pnext"].append(x.next_reference_start)
                        df_dict["roi"].append(row["roi"])
                except:
                    pass
            end_list = np.core.defchararray.add(
                "{}~{}__".format(row["roi"], str(row["end_pos"])),
                np.arange(row["start_count"], row["start_count"] + row["end_count"] + 1, 1)
                .astype(int)
                .astype(str),
            )
            for read in end_list:
                try:
                    end_reads = names.find(read)
                    for x in end_reads:
                        df_dict["qname"].append(x.query_name)
                        df_dict["rname"].append(x.reference_name)
                        df_dict["pos"].append(x.reference_start)
                        df_dict["cigar"].append(x.cigarstring)
                        df_dict["tlen"].append(x.template_length)
                        df_dict["seq"].append(x.query_sequence)
                        df_dict["flag"].append(x.flag)
                        df_dict["pnext"].append(x.next_reference_start)
                        df_dict["roi"].append(row["roi"])
                except:
                    pass
        df = pd.DataFrame.from_dict(df_dict)
        df = df.drop_duplicates(
            subset=["qname", "rname", "pos", "cigar", "tlen", "seq", "flag", "pnext"],
            keep="first",
        )
        df = df.dropna(subset=["rname"])
        df["tlen"] = df["tlen"].astype(int)
        df["pos"] = df["pos"].astype(int)
        self.clip_end_df = df


    def get_clip_pos_end_rois(
        self,
        raw_seq_dict: Dict[str, Any]
    ) -> None:
        """
        Get positions for clipped ends (end rois).

        Args:
            raw_seq_dict (Dict[str, Any]): Dictionary of sequences. 

        Returns:
            None
        """

        self.clip_end_df = self.clip_end_df[
            (self.clip_end_df["cigar"].str.count("M") == 1)
            & (self.clip_end_df["cigar"].str.count("S") == 0)
            & (self.clip_end_df["cigar"].str.count("I") == 0)
            & (self.clip_end_df["cigar"].str.count("X") == 0)
            & (self.clip_end_df["cigar"].str.count("D") == 0)
        ]

        sub_df = (
            self.clip_end_df["qname"]
            .str.split("__")
            .str[0]
            .value_counts()
            .rename_axis("coords")
            .reset_index(name="counts")
        )

        df_list = []

        for index, row in self.clip_end_df.iterrows():
            self.clip_end_df.loc[index, "slen"] = len(row["seq"])

        for index, row in self.end_roi_df.iterrows():
            sub_list = []
            start_name = row["roi"] + "~" + str(row["start_pos"])
            end_name = row["roi"] + "~" + str(row["end_pos"])

            start_df = self.clip_end_df[self.clip_end_df["qname"].str.split("__").str[0] == start_name].copy()
            end_df = self.clip_end_df[self.clip_end_df["qname"].str.split("__").str[0] == end_name].copy()

            start_df["query"] = start_df["qname"].str.split("__").str[0].copy()
            start_df = start_df.groupby(["query", "rname"])["qname"].count()
            start_df = start_df[start_df >= 10]

            end_df["query"] = end_df["qname"].str.split("__").str[0].copy()
            end_df = end_df.groupby(["query", "rname"])["qname"].count()
            end_df = end_df[end_df >= 10]

            if (len(start_df) > 0) and (len(end_df) > 0):
                start_df = start_df.reset_index(name="counts")

                for i in start_df["rname"].unique():
                    sub_df = self.clip_end_df[
                        (self.clip_end_df["qname"].str.split("__").str[0] == start_name)
                        & (self.clip_end_df["rname"] == i)
                    ].copy()
                    start_contig_length = len(raw_seq_dict[i].seq)
                    sub_df2 = sub_df.groupby(
                        pd.cut(
                            sub_df["pos"],
                            np.arange(-1, start_contig_length + 201, 200)
                        )
                    )["qname"].count()
                    sub_df2 = sub_df2[sub_df2 >= 10]
                    sub_df2 = sub_df2.reset_index(name="counts")
                    sub_df2["start_query_pos"] = (
                        start_df["query"].str.split("~").str[-1].astype(int)
                    )
                    sub_df2["query"] = (
                        start_df["query"].str.split("~").str[:4].str.join("~")
                    )
                    sub_df2 = sub_df2.sort_values(by=["counts"], ascending=False).iloc[[0]]

                    for ind, r in sub_df2.iterrows():
                        if not isinstance(r["pos"], float):
                            sub_df2.loc[ind, "left"] = r["pos"].left
                            sub_df2.loc[ind, "right"] = r["pos"].right
                        else:
                            sub_df2.loc[ind, "left"] = np.nan
                            sub_df2.loc[ind, "right"] = np.nan

                    sub_df = sub_df[
                        sub_df["pos"].between(
                            sub_df2["left"].iloc[0], sub_df2["right"].iloc[0]
                        )
                    ]
                    self.end_roi_df.loc[index, "right_end_pos"] = np.median(sub_df["pos"])
                    self.end_roi_df.loc[index, "right_end_count"] = sub_df2["counts"].iloc[0]
                    self.end_roi_df.loc[index, "right_end_rname"] = i

                end_df = end_df.reset_index(name="counts")

                for i in end_df["rname"].unique():
                    sub_df = self.clip_end_df[
                        (self.clip_end_df["qname"].str.split("__").str[0] == end_name)
                        & (self.clip_end_df["rname"] == i)
                    ].copy()
                    end_contig_length = len(raw_seq_dict[i].seq)
                    sub_df2 = sub_df.groupby(
                        pd.cut(
                            sub_df["pos"],
                            np.arange(-1, end_contig_length + 201, 200)
                        )
                    )["qname"].count()
                    sub_df2 = sub_df2[sub_df2 >= 10]
                    sub_df2 = sub_df2.reset_index(name="counts")
                    sub_df2["end_query_pos"] = (
                        end_df["query"].str.split("~").str[-1].astype(int)
                    )
                    sub_df2["query"] = end_df["query"].str.split("~").str[:4].str.join("~")
                    sub_df2 = sub_df2.sort_values(by=["counts"], ascending=False).iloc[[0]]

                    for ind, r in sub_df2.iterrows():
                        sub_df2.loc[ind, "left"] = r["pos"].left
                        sub_df2.loc[ind, "right"] = r["pos"].right

                    sub_df = sub_df[
                        sub_df["pos"].between(
                            sub_df2["left"].iloc[0], sub_df2["right"].iloc[0]
                        )
                    ]
                    self.end_roi_df.loc[index, "left_end_pos"] = np.median(sub_df["pos"])
                    self.end_roi_df.loc[index, "left_end_count"] = sub_df2["counts"].iloc[0]
                    self.end_roi_df.loc[index, "left_end_rname"] = i

            else:
                self.end_roi_df.loc[index, "right_end_pos"] = np.nan
                self.end_roi_df.loc[index, "right_end_count"] = np.nan
                self.end_roi_df.loc[index, "right_end_rname"] = np.nan
                self.end_roi_df.loc[index, "left_end_pos"] = np.nan
                self.end_roi_df.loc[index, "left_end_count"] = np.nan
                self.end_roi_df.loc[index, "left_end_rname"] = np.nan

        self.end_roi_df = self.end_roi_df.dropna()

        if len(self.end_roi_df) < 1:
            logger.into("No ROIs spanning the ends with clipped reads.")



    def filter_by_end_status(self) -> pd.DataFrame:
        """
        Filter end rois based on end_info and contig matching.

        Returns:
            pd.DataFrame: DataFrame with filtered ROI information.
        """

        
        for index, row in self.end_roi_df.iterrows():
            if (row["end_info"] == "whole") and (
                (row["contig"] == row["right_end_rname"])
                or (row["contig"] == row["left_end_rname"])
            ):
                self.end_roi_df.drop(index, inplace=True)
            elif (row["end_info"] == "left") and (row["contig"] == row["left_end_rname"]):
                self.end_roi_df.drop(index, inplace=True)
            elif (row["end_info"] == "right") and (row["contig"] == row["right_end_rname"]):
                self.end_roi_df.drop(index, inplace=True)
        
        if len(self.end_roi_df) < 1:
            logger.info("No ROIs spanning the ends with clipped reads.")

        


    def reformat_end_roi_tables(self, seq: Dict[str, Any]) -> None:
        """
        Reformat end ROI tables based on end_info and sequence information.

        Args:
            seq (Dict[str, Any]): A dictionary containing sequence information.

        Returns:
            None
        """

        (
            self.end_roi_df["start_pos"],
            self.end_roi_df["end_pos"],
            self.end_roi_df["left_end_pos"],
            self.end_roi_df["right_end_pos"],
        ) = (
            self.end_roi_df["start_pos"].astype(int),
            self.end_roi_df["end_pos"].astype(int),
            self.end_roi_df["left_end_pos"].astype(int),
            self.end_roi_df["right_end_pos"].astype(int),
        )
        for index, row in self.end_roi_df.iterrows():
            if row["end_info"] == "whole":
                self.end_roi_df = fix_table_whole(index, row, self.end_roi_df, seq)
            elif row["end_info"] == "left":
                self.end_roi_df = fix_table_left(index, row, self.end_roi_df, seq)
            elif row["end_info"] == "right":
                self.end_roi_df = fix_table_right(index, row, self.end_roi_df, seq)
        self.end_roi_df = self.end_roi_df[
            [
                "roi",
                "start_pos",
                "start_count",
                "end_pos",
                "end_count",
                "length",
                "contig",
                "contig_len",
                "end_info",
                "contig_split",
                "med_z",
            ]
        ]

    """
    long
    """


    def bamfile_to_pandas(self, output: Path, neighbourhood: int) -> None:
        """
        Parse and process a sorted BAM file into a Pandas DataFrame.

        Parameters:
            output (Path): The output folder where the SAM file is located.
            neighbourhood (int): number of bases before and after the ROI
        
        Returns:
            None

        This method parses and processes a sorted BAM file into a Pandas DataFrame. It filters reads based on the ROIs provided
        in 'self.roi_df', ensuring that only reads within the specified regions are included in the output DataFrame.
        
        """

        sorted_bam: Path = output / "temp_minimap_sorted.bam"

        bamfile = pysam.AlignmentFile(sorted_bam, "rb")

        df_dict: Dict[str, list] = {
            "qname": [],
            "rname": [],
            "pos": [],
            "cigar": [],
            "tlen": [],
            "seq": [],
            "flag": [],
            "pnext": [],
            "qqual": [],
        }

        for index, row in self.roi_df.iterrows():
            contig = row["accession"].split("~")[0]
            start = row["start"] - neighbourhood
            if start < 0:
                start = 0
            end = row["end"] + neighbourhood
            if end > row["contig_len"]:
                end = row["contig_len"]
            for reads in bamfile.fetch(contig, start, end):
                df_dict["qname"].append(reads.query_name)
                df_dict["rname"].append(contig)
                df_dict["pos"].append(reads.reference_start)
                df_dict["cigar"].append(reads.cigarstring)
                df_dict["tlen"].append(reads.template_length)
                df_dict["seq"].append(reads.query_sequence)
                df_dict["flag"].append(reads.flag)
                df_dict["pnext"].append(reads.next_reference_start)
                df_dict["qqual"].append(reads.query_qualities.tolist())

        df = pd.DataFrame.from_dict(df_dict)
        df = df.drop_duplicates(
            subset=["qname", "rname", "pos", "cigar", "tlen", "seq", "flag", "pnext"],
            keep="first",
        )
        df["tlen"] = df["tlen"].astype(int)
        df["pos"] = df["pos"].astype(int)

        if df.empty:
            exit_error_gracefully(self.output, self.all_zscores, self.depths, self.median, self.mad, self.start_time)
            logger.error("No depth found at all positions.")

        self.sam_df = df

#### do i need the secondaries?

    def bamfile_to_pandas_long(self, output: Path, neighbourhood: int) -> None:
        """
        Parse and process a sorted BAM file into a Pandas DataFrame for long reads

        Parameters:
            output (Path): The output folder where the SAM file is located.
            neighbourhood (int): number of bases before and after the ROI
        
        Returns:
            None

        This method parses and processes a sorted BAM file into a Pandas DataFrame. It filters reads based on the ROIs provided
        in 'self.roi_df', ensuring that only reads within the specified regions are included in the output DataFrame.
        """

        sorted_bam: Path = output / "temp_minimap_sorted.bam"

        bamfile = pysam.AlignmentFile(sorted_bam, "rb")

        df_dict: Dict[str, list] = {
            'qname':[], 'rname':[], 'pos':[], 'cigar':[],'tlen':[],'seq':[],'flag':[],'pnext':[], 'qqual':[], 'reference_end':[], 'reference_start': []
        }
        df_dict_secondary: Dict[str, list] = {
            'qname':[], 'rname':[], 'pos':[], 'cigar':[],'tlen':[],'seq':[],'flag':[],'pnext':[], 'qqual':[], 'reference_end':[], 'reference_start': []
        }

        primary_flags = [0, 16, 2048, 2064]
        secondary_flags = [256, 272]
        for index, row in self.roi_df.iterrows():
            contig = row['accession'].split('~')[0]
            start = row['start'] - neighbourhood
            if start < 0:
                start = 0
            end = row['end'] + neighbourhood
            if end > row['contig_len']:
                end = row['contig_len']
        # ONT reads - issue with non-primary alignments
        # keep only 0 (mapped fwd) or 16 (mapped rev)
            for reads in bamfile.fetch(contig, start, end):
                if reads.flag in primary_flags: # keeps all the primary and supp reads
                    df_dict['qname'].append(reads.query_name)
                    df_dict['rname'].append(contig)
                    df_dict['pos'].append(reads.reference_start)
                    df_dict['cigar'].append(reads.cigarstring)
                    df_dict['tlen'].append(reads.template_length)
                    df_dict['seq'].append(reads.query_sequence)
                    df_dict['flag'].append(reads.flag)
                    df_dict['pnext'].append(reads.next_reference_start)
                    df_dict['reference_end'].append(reads.reference_end)
                    df_dict['reference_start'].append(reads.reference_start)
                    df_dict['qqual'].append(reads.query_qualities.tolist())
                elif reads.flag in primary_flags:
                    df_dict_secondary['qname'].append(reads.query_name)
                    df_dict_secondary['rname'].append(contig)
                    df_dict_secondary['pos'].append(reads.reference_start)
                    df_dict_secondary['cigar'].append(reads.cigarstring)
                    df_dict_secondary['tlen'].append(reads.template_length)
                    df_dict_secondary['seq'].append(reads.query_sequence)
                    df_dict_secondary['flag'].append(reads.flag)
                    df_dict_secondary['pnext'].append(reads.next_reference_start)
                    df_dict_secondary['reference_end'].append(reads.reference_end)
                    df_dict_secondary['reference_start'].append(reads.reference_start)
                    df_dict_secondary['qqual'].append(reads.query_qualities.tolist())

        df = pd.DataFrame.from_dict(df_dict)
        df = df.drop_duplicates(subset = ['qname', 'rname', 'pos', 'cigar', 'tlen', 'seq', 'flag', 'pnext', 'reference_end', 'reference_start'], keep='first')
        df['tlen'] = df['tlen'].astype(int)
        df['pos'] = df['pos'].astype(int)

        sam_secondary_df = pd.DataFrame.from_dict(df_dict_secondary)
        sam_secondary_df = sam_secondary_df.drop_duplicates(subset = ['qname', 'rname', 'pos', 'cigar', 'tlen', 'seq', 'flag', 'pnext', 'reference_end', 'reference_start'], keep='first')
        sam_secondary_df['tlen'] = sam_secondary_df['tlen'].astype(int)
        sam_secondary_df['pos'] = sam_secondary_df['pos'].astype(int)

        if df.empty:
            exit_error_gracefully(self.output, self.all_zscores, self.depths, self.median, self.mad, self.start_time)
            logger.error("No depth found at all positions.")
        # not big a deal if secondary df empty
        # if sam_secondary_df.empty: 
        #     logger.warning("No secondary alignment depth found at all positions.")

        self.sam_df = df
        self.sam_secondary_df = sam_secondary_df

 
    def find_multimapped_reads_long(self, neighbourhood: int, min_reads: int) -> None:
        """
        Identify ROIs with multiple mapped long reads that span both ends.

        Args:
            neighbourhood (int): The size of the neighbourhood around ROI start and end positions.
            min_reads (int): The minimum number of long reads required to consider an ROI.

        Returns:
            None
        """

        # gets all the reads that have more than 1 primary or supplementary mapping
        df = self.sam_df[self.sam_df.duplicated('qname', keep=False) == True]
        # loop over each roi
        for index, row in self.roi_df.iterrows():
            center = row['center']
            # gets all the reads in the vicinity of the start and end of the ROI
            sub_df = df[(df['rname'] == row['contig']) &
                        (((df['pos'].between(row['start'] - neighbourhood, center)) | (df['pnext'].between(center, row['end'] + neighbourhood))) |
                        ((df['pnext'].between(row['start'] - neighbourhood, center)) | (df['pos'].between(center, row['end'] + neighbourhood))))].copy()
            # keep only the dupes again 
            sub_df = sub_df[sub_df.duplicated('qname', keep=False) == True]
            # collect all starts and ends
            # automatically sorts the list from largest to smalles
            starts = sub_df['reference_start'].value_counts(dropna=False)
            ends = sub_df['reference_end'].value_counts(dropna=False)
            # Don't keep empty ones
            if (len(sub_df)/2) < 1:
                self.roi_df = self.roi_df[self.roi_df['accession'] != row['accession']].copy()
            else:
                # get the top hit
                start_coverage = int(starts.iloc[0])
                start_position = int(starts.index[0])
                end_coverage = int(ends.iloc[0])
                end_position = int(ends.index[0])
                # need more than min_reads coverage at the start and end (default min_reads = 5)
                # skips the row if there is less than 5 at either end
                if start_coverage < min_reads or end_coverage < min_reads:
                    self.roi_df = self.roi_df[self.roi_df['accession'] != row['accession']].copy()
                else:
                # denotes the read with likely_start and likely_end
                    self.roi_df.loc[index,'likely_start'] = start_position
                    self.roi_df.loc[index,'likely_end'] = end_position
                    self.roi_df.loc[index,'start_count'] = start_coverage
                    self.roi_df.loc[index,'end_count'] = end_coverage
        if len(self.roi_df) < 1:
            exit_error_gracefully(self.output, self.all_zscores, self.depths, self.median, self.mad, self.start_time)
            logger.exit(f"No ROIs found with at least {min_reads} long reads that span both ends.")
        else:
            logger.info(f"{len(self.roi_df)} ROI(s) found with at least {min_reads} long reads that span both ends.")
        # return the roi_df







"""


Auxillary functions (from original hafez)



"""


def check_read_lengths(row, index, clip_df, side_name, names, side):
    df_dict = {
        "qname": [],
        "rname": [],
        "pos": [],
        "cigar": [],
        "tlen": [],
        "seq": [],
        "flag": [],
        "pnext": [],
        "roi": [],
        "seq_len": [],
    }
    
    start_list = np.core.defchararray.add(
        "{}~{}__".format(row["roi"], str(row["start_pos"])),
        np.arange(1, row["start_count"] + 1, 1).astype(int).astype(str),
    )
    if side == "start":
        for read in start_list:
            try:
                start_reads = names.find(read)
                for x in start_reads:
                    df_dict["qname"].append(x.query_name)
                    df_dict["rname"].append(x.reference_name)
                    df_dict["pos"].append(x.reference_start)
                    df_dict["cigar"].append(x.cigarstring)
                    df_dict["tlen"].append(x.template_length)
                    df_dict["seq"].append(x.query_sequence)
                    df_dict["flag"].append(x.flag)
                    df_dict["pnext"].append(x.next_reference_start)
                    df_dict["roi"].append(row["roi"])
                    df_dict["seq_len"].append(len(x.query_sequence))
            except:
                pass
        start_df = pd.DataFrame.from_dict(df_dict)
        start_df = start_df.drop_duplicates()
        if (
            len(start_df[start_df["seq_len"] >= 35]) < 10
        ):  ##### this is to make sure that any reads that didnt have clipped sections that were no long enough to be mapped are not unfairly excluded from the search
            return 10
    elif side == "end":
        end_list = np.core.defchararray.add(
            "{}~{}__".format(row["roi"], str(row["end_pos"])),
            np.arange(row["start_count"], row["start_count"] + row["end_count"] + 1, 1)
            .astype(int)
            .astype(str),
        )
        for read in end_list:
            try:
                end_reads = names.find(read)
                for x in end_reads:
                    df_dict["qname"].append(x.query_name)
                    df_dict["rname"].append(x.reference_name)
                    df_dict["pos"].append(x.reference_start)
                    df_dict["cigar"].append(x.cigarstring)
                    df_dict["tlen"].append(x.template_length)
                    df_dict["seq"].append(x.query_sequence)
                    df_dict["flag"].append(x.flag)
                    df_dict["pnext"].append(x.next_reference_start)
                    df_dict["roi"].append(row["roi"])
                    df_dict["seq_len"].append(len(x.query_sequence))
            except:
                pass
        end_df = pd.DataFrame.from_dict(df_dict)
        end_df = end_df.drop_duplicates()
        if (
            len(end_df[end_df["seq_len"] >= 35]) < 10
        ):  ##### this is to make sure that any reads that didnt have clipped sections that were no long enough to be mapped are not unfairly excluded from the search
            return 10



def fix_table_whole(index, row, roi_df, seq):
    left_contig_length = len(seq[row["left_end_rname"]].seq)
    closest_left = min(
        [0, left_contig_length], key=lambda x: abs(x - row["left_end_pos"])
    )
    if closest_left == 0:
        left_start = row["left_end_pos"]
        left_end = 0
    else:
        left_start = left_contig_length
        left_end = row["left_end_pos"]
    right_contig_length = len(seq[row["right_end_rname"]].seq)
    closest_right = min(
        [0, right_contig_length], key=lambda x: abs(x - row["right_end_pos"])
    )
    if closest_right == 0:
        right_start = 0
        right_end = row["right_end_pos"]
    else:
        right_start = row["right_end_pos"]
        right_end = righ_contig_length
    roi_df.loc[index, "contig_split"] = "{}({}..{}) -> {}({}..{}) -> {}({}..{})".format(
        row["left_end_rname"],
        left_start,
        left_end,
        row["contig"],
        row["start_pos"],
        row["end_pos"],
        row["right_end_rname"],
        right_start,
        right_end,
    )
    return roi_df


def fix_table_left(index, row, roi_df, seq):
    left_contig_length = len(seq[row["left_end_rname"]].seq)
    closest_left = min(
        [0, left_contig_length], key=lambda x: abs(x - row["left_end_pos"])
    )
    if closest_left == 0:
        left_start = row["left_end_pos"]
        left_end = 0
    else:
        left_start = left_contig_length
        left_end = row["left_end_pos"]
    roi_df.loc[index, "contig_split"] = "{}({}..{}) -> {}({}..{})".format(
        row["left_end_rname"],
        left_start,
        left_end,
        row["contig"],
        row["start_pos"],
        row["end_pos"],
    )
    return roi_df


def fix_table_right(index, row, roi_df, seq):
    right_contig_length = len(seq[row["right_end_rname"]].seq)
    closest_right = min(
        [0, right_contig_length], key=lambda x: abs(x - row["right_end_pos"])
    )
    if closest_right == 0:
        right_start = row["right_end_pos"]
        right_end = 0
    else:
        right_start = right_contig_length
        right_end = row["right_end_pos"]
    roi_df.loc[index, "contig_split"] = "{}({}..{}) -> {}({}..{})".format(
        row["contig"],
        row["start_pos"],
        row["end_pos"],
        row["right_end_rname"],
        right_start,
        right_end,
    )
    return roi_df



def get_longest_len_z(z, cutoff):
    x = {"start": [], "end": []}
    above_threshold = np.diff(
        np.array(z) < cutoff, prepend=False
    )  # less than as this then will flag start and end regions as where the go below and come above
    rise_above_thres = np.argwhere(above_threshold)[::2, 0]
    sink_below_thres = np.argwhere(above_threshold)[1::2, 0]
    x["start"] = rise_above_thres
    x["end"] = sink_below_thres
    lengths = []
    if len(x["start"]) > 0:
        for j in np.arange(0, len(x["start"]) - 1):
            lengths.append(x["end"][j] - x["start"][j])
        if len(lengths) > 0:
            longest_below_z = max(lengths)
        else:
            longest_below_z = 0
    else:
        longest_below_z = 0
    return longest_below_z


def multiprocessing_calc_roi_z_medians(roi_df, depths, median, mad, median_z_cutoff):
    for i in roi_df["contig"].unique():
        z = [((0.6745 * (x - median)) / mad) for x in depths[i]]
        for index, row in roi_df[(roi_df["contig"] == i)].iterrows():
            med_z = z[row["start_pos"] : row["end_pos"]]
            med_z = np.median(med_z)
            med_z_l = z[(row["start_pos"] - 10000) : row["start_pos"]]
            longest_below_z = get_longest_len_z(
                z[row["start_pos"] : row["end_pos"]], (med_z / 2)
            )
            if len(med_z_l) > 0:
                med_z_l = np.median(med_z_l)
            else:
                med_z_l = np.nan
            med_z_r = z[row["end_pos"] : (row["end_pos"] + 10000)]
            if len(med_z_r) > 0:
                med_z_r = np.median(med_z_r)
            else:
                med_z_r = np.nan
            roi_df.loc[index, "med_z"] = med_z
            roi_df.loc[index, "med_z_l"] = med_z_l
            roi_df.loc[index, "med_z_r"] = med_z_r
            roi_df.loc[index, "longest_below_z"] = longest_below_z
    return roi_df









def multiprocess_find_soft_clippings(roi_df: pd.DataFrame, clips: pd.DataFrame, width: int, neighbourhood: int) -> pd.DataFrame:
    """
    Find soft clippings within ROIs and create a DataFrame with the results.

    For find_soft_clippings multiprocess parsing

    Parameters:
        width (int): The minimum width threshold for soft clippings to be considered.
        neighbourhood (int): number of bases before and after the ROI

    Returns:
       df

    This method identifies soft clippings within the ROIs in 'roi_df' and creates a DataFrame
    with information about their positions and counts. It filters the DataFrame based on the provided width
    threshold and returns it.

    """

    clip_dict = {
        "roi": [],
        "start_pos": [],
        "start_count": [],
        "end_pos": [],
        "end_count": [],
        "length": [],
        "contig_len": [],
        "contig": [],
        "end_info": [],
    }
    for index, row in roi_df.iterrows():
        s_list = []
        sub_df = clips[
            (clips["rname"] == row["contig"])
            & (clips["pos"].between(row["start"] - neighbourhood, row["end"] + neighbourhood))
        ]
        for i, r in sub_df.iterrows():
            M = re.sub("[0-9]+", "", r["cigar"]).find("M")
            S = re.sub("[0-9]+", "", r["cigar"]).find("S")
            if S == 1:
                s_list.append(r["pos"] + int(re.split(r"\D+", r["cigar"])[M]))
            else:
                s_list.append(r["pos"])
        if len(s_list) > 1:
            sub_df = (
                pd.Series(s_list)
                .value_counts()
                .rename_axis("coords")
                .reset_index(name="counts")
            )
            sub_df = sub_df[sub_df["counts"] >= 10]
            if (len(sub_df[sub_df["coords"] < row["center"]]) > 0) and (
                len(sub_df[sub_df["coords"] > row["center"]]) > 0
            ):
                for i, r in sub_df[sub_df["coords"] < row["center"]].iterrows():
                    for ind, rw in sub_df[sub_df["coords"] > row["center"]].iterrows():
                        clip_dict["roi"].append(row["accession"])
                        clip_dict["start_pos"].append(r["coords"])
                        clip_dict["start_count"].append(r["counts"])
                        clip_dict["end_pos"].append(rw["coords"])
                        clip_dict["end_count"].append(rw["counts"])
                        clip_dict["length"].append(rw["coords"] - r["coords"])
                        clip_dict["contig_len"].append(row["contig_len"])
                        clip_dict["contig"].append(row["contig"])
                        clip_dict["end_info"].append(row["end_info"])
    df = pd.DataFrame.from_dict(clip_dict)
    df = df[df["length"] >= width]
    return df


def collect_clips(roi_df, clips):
    read_list = []
    for index, row in roi_df.iterrows():
        sub_df = clips[
            (clips["rname"] == row["contig"])
            & (clips["pos"].between(row["start_pos"] - 200, row["end_pos"] + 200))
        ]
        counter_2 = 1
        for i, r in sub_df.iterrows():
            M = re.sub("[0-9]+", "", r["cigar"]).find("M")
            S = re.sub("[0-9]+", "", r["cigar"]).find("S")
            if S == 1:
                matched = int(re.split(r"\D+", r["cigar"])[M])
                pos = r["pos"] + matched
                clipped_seq = r["seq"][matched:]
                clipped_qual = r["qqual"][matched:]
            else:
                matched = int(re.split(r"\D+", r["cigar"])[S])
                pos = r["pos"]
                clipped_seq = r["seq"][:matched]
                clipped_qual = r["qqual"][:matched]
            if (pos == row["start_pos"]) or (pos == row["end_pos"]):
                clip_name = row["roi"] + "~" + str(pos) + "__" + str(counter_2)
                clipped = SeqRecord(
                    Seq(clipped_seq),
                    id=clip_name,
                    name=clip_name,
                    description="",
                    dbxrefs=[],
                )
                clipped.letter_annotations["phred_quality"] = clipped_qual
                read_list.append(clipped)
                counter_2 = counter_2 + 1
    return [read_list, roi_df]
