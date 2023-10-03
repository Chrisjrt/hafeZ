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
from hafeZ.utils.post_processing import output_contig_Z
from hafeZ.utils.util import remove_file


def exit_error_gracefully(
    output: Path,
    all_zscores: bool,
    depths: Dict[Union[str, int], List[float]],
    median: float,
    mad: float,
    start_time: float,
) -> None:
    """
    Gracefully handle errors and perform cleanup tasks.

    Args:
        output (Path): The output directory.
        all_zscores (bool): Flag indicating whether to calculate and save z-scores.
        depths (dict): A dictionary of depth data.
        median (float): The median value.
        mad (float): The Median Absolute Deviation (MAD) value.
        start_time (float): The start time of the process.

    Returns:
        None
    """
    # output empty roi summary
    output_no_roi(output)
    # remove temp files
    clean_up(output)
    if all_zscores is True:
        output_contig_Z(depths, output, median, mad)

    # Determine elapsed time
    elapsed_time = time.time() - start_time
    elapsed_time = round(elapsed_time, 2)

    # Show elapsed time for the process
    logger.info("hafeZ run complete!")
    logger.info(f"Total run time: {elapsed_time} seconds.")


def exit_error_gracefully_premap(output: Path, start_time: float) -> None:
    """
    Gracefully handle errors and perform cleanup tasks in Premap class.

    Args:
        output (Path): The output directory.
        start_time (float): The start time of the process.

    Returns:
        None
    """
    # output empty roi summary
    output_no_roi(output)
    # remove temp files
    clean_up(output)

    # Determine elapsed time
    elapsed_time = time.time() - start_time
    elapsed_time = round(elapsed_time, 2)

    # Show elapsed time for the process
    logger.info("hafeZ run complete!")
    logger.info(f"Total run time: {elapsed_time} seconds.")


def exit_success(output: Path, start_time: float) -> None:
    """
    Handle successful exit and perform cleanup tasks.

    Args:
        output (Path): The output directory.
        start_time (float): The start time of the process.

    Returns:
        None
    """

    clean_up(output)

    # Determine elapsed time
    elapsed_time = time.time() - start_time
    elapsed_time = round(elapsed_time, 2)

    # Show elapsed time for the process
    logger.info("hafeZ run complete!")
    logger.info(f"Total run time: {elapsed_time} seconds.")


def output_no_roi(output: Path) -> None:
    """
    Output an empty ROI summary and log a message about no active prophage found.

    Args:
        output (Path): The output directory.

    Returns:
        None
    """
    roi_df = pd.DataFrame(
        data={
            "roi": [np.nan],
            "start_pos": [np.nan],
            "end_pos": [np.nan],
            "roi_length": [np.nan],
            "orf_count": [np.nan],
            "frac_phrog": [np.nan],
            "circular": [np.nan],
        }
    )
    out_file: Path = Path(output) / "hafeZ_summary_all_rois.tsv"
    roi_df.to_csv(out_file, sep="\t", index=False)
    logger.info("NO ACTIVE PROPHAGES FOUND.")
    logger.info(f"Exiting and outputing an empty .tsv summary file {out_file}.")


def clean_up(output: Path) -> None:
    """
    Clean up temporary files in the specified output directory.

    Args:
        output (Path): The output directory.

    Returns:
        None
    """
    for temp_file in glob.glob(str(output) + "/temp*"):
        remove_file(temp_file)
