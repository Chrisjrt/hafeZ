"""
includes all process_rois.py functions from v1 that aren't in the haf class
"""

import itertools
import multiprocessing
import os
import re
import subprocess
import sys
import time
from collections import Counter
from functools import partial
from io import StringIO
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from loguru import logger
from scipy import stats

from hafez.utils.external_tools import ExternalTool
from hafez.utils.util import split_df_for_multiproc


def quick_filter(roi_df: pd.DataFrame) -> pd.DataFrame:
    """
    Quickly filter rows in the ROI DataFrame based on specified criteria.

    Args:
        roi_df (pd.DataFrame): The DataFrame containing ROI information.

    Returns:
        pd.DataFrame: Filtered DataFrame with selected rows.
    """
    df_list = []
    for i in roi_df["roi"].unique():
        sub_df = roi_df[roi_df["roi"] == i].copy()
        for index, row in sub_df.iterrows():
            if row["circular"] != True:
                sub_df.loc[index, "total_clips"] = int(row["start_count"]) + int(
                    row["end_count"]
                )
            else:
                sub_df.loc[index, "total_clips"] = None  # Use None for NaN values
        sub_df = sub_df.sort_values(by=["total_clips"], ascending=False).iloc[[0]]
        df_list.append(sub_df)
    roi_df = pd.concat(df_list)
    roi_df = roi_df.reset_index(drop=True)
    roi_df["contig_split"] = None  # Use None for NaN values

    return roi_df
