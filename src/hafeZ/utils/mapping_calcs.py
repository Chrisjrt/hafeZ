import time
import pandas as pd
from scipy.signal import savgol_filter
from scipy.stats import zscore
from scipy.signal import find_peaks
import numpy as np
from pathlib import Path
import sys
import matplotlib.pyplot as plt
from statistics import mean
from statistics import stdev
import multiprocessing
from functools import partial
from pandas import DataFrame
from typing import Dict, List, Union, Tuple
from loguru import logger
from hafeZ.utils.exit import exit_error_gracefully_premap

def savgol(df: DataFrame, bin_size: int) -> List[Dict[Union[str, int], List[float]]]:
    """
    Apply a Savitzky-Golay filter to smooth depth data in a DataFrame.

    Parameters:
        df (DataFrame): A DataFrame containing depth data, typically with columns 'acc' and 'depth'.
        bin_size (int): The window size for the Savitzky-Golay filter.

    Returns:
        List[Dict[Union[str, int], List[float]]]: A list of dictionaries containing raw and smoothed depth data.

    This function takes a DataFrame 'df' with depth data and applies a Savitzky-Golay filter to smooth the data.
    It calculates both raw and smoothed depth data for each unique 'acc' value in the DataFrame. The results are
    returned as a list of dictionaries, where each dictionary has keys 'acc' and 'depth' representing the
    account identifier and the corresponding depth values.

    The 'bin_size' parameter controls the size of the moving window used by the Savitzky-Golay filter.

    Example:
        df = pd.DataFrame({'acc': ['A', 'A', 'B', 'B'], 'depth': [10, 12, 8, 7]})
        smoothed_data = savgol(df, 3)
    """
    df = df[1]
    raw_dict = {}
    smoothed_dict = {}
    raw_dict[df["acc"].unique()[0]] = df["depth"]
    z = savgol_filter(df["depth"], window_length=bin_size,  polyorder=2)
    smoothed_dict[df["acc"].unique()[0]] = z
    return [raw_dict, smoothed_dict]




def smooth_depths(output: Path, bin_size: int, threads: int) -> Tuple[Dict[Union[str, int], List[float]], Dict[Union[str, int], List[float]]]:
    """
    Smooth depth data from BED files using a Savitzky-Golay filter.

    Parameters:
        output_folder (str): Path to the output directory where the BED file is located.
        bin_size (int): The window size for the Savitzky-Golay filter.
        threads (int): Number of threads/cores to use for parallel processing.

    Returns:
        Tuple[Dict[Union[str, int], List[float]], Dict[Union[str, int], List[float]]]: A tuple containing dictionaries
        with smoothed and raw depth data.

    This function reads a BED file containing depth data from the specified 'output', groups the data by 'acc' values,
    and applies a Savitzky-Golay filter to smooth the depth data. The smoothed and raw depth data are returned as two dictionaries.

    The 'bin_size' parameter controls the size of the moving window used by the Savitzky-Golay filter, and 'threads' specifies
    the number of CPU cores for parallel processing.

    Example:
        smoothed, raw = smooth_depths("output_dir", 5, 4)
    """
    
    # Read the BED file into a DataFrame

    bed_file: Path = output / "temp_minimap.regions.bed.gz"

    df = pd.read_csv(
        bed_file,
        compression="gzip",
        sep="\t",
        names=("acc", "start", "end", "depth"),
    )


    # Group data by 'acc' values
    df = df.groupby("acc")

    # Create a multiprocessing pool for parallel processing
    pool = multiprocessing.Pool(processes=threads)

    ####
    # add checks for bin_size > smallest contig here
    ########

    # Apply the Savitzky-Golay filter to each group and store the results in a list of dictionaries
    dicts = pool.map(partial(savgol, bin_size=bin_size), df)

    # Merge dictionaries into a single dictionary for raw and smoothed data
    raw_dict = {}
    smoothed_dict = {}
    for i in dicts:
        acc = list(i[0].keys())[0]
        raw_dict[acc] = list(i[0].values())[0]
        smoothed_dict[acc] = list(i[1].values())[0]

    return smoothed_dict, raw_dict


def plot_MAD_error_coverage(cov: Dict[Union[str, int], List[float]], output: Path) -> None:
    """
    Plot coverage depth for each contig and save the plots to the output directory.

    Args:
        cov (Dict[Union[str, int], List[float]]): A dictionary where keys are contig names or IDs,
            and values are lists of coverage depth values.
        output (Path): The path to the output directory where the plots will be saved.

    Returns:
        None
    """
    for i in cov:
        depth_to_plot = cov[i]
        fig = plt.gcf()
        fig.set_size_inches(18.5, 10.5)
        plt.plot(depth_to_plot)
        plt.ylabel('Coverage depth')
        plt.savefig(output / ('MAD_ERROR_coverage_plot_for_' + str(i) + '.png'), format='png')
        plt.clf()
    print('\n{:#^50}'.format(''))



def get_ZScores(depths: Dict[Union[str, int], List[float]], output: Path, start_time: float, cov: Dict[Union[str, int], List[float]], expect_mad_zero: bool) -> Tuple[Dict[Union[str, int], List[float]], float, float]:
    """
    Calculate Z-scores for depth data.

    Parameters:
        depths (Dict[Union[str, int], List[float]]): A dictionary containing depth data for each contig or region.
        output (Path): output directory.
        start_time (float): start time.
        cov (Dict[Union[str, int], List[float]]): A dictionary where keys are contig names or IDs,
            and values are lists of coverage depth values.
        expect_mad_zero (bool): Will also cause coverage plots for each contig to be output to help with debugging. Useful for uninduced lysates

    Returns:
        Tuple[Dict[Union[str, int], List[float]], float, float]: A tuple containing Z-scores for each contig/region,
        the median of the depth data, and the Median Absolute Deviation (MAD).

    This function calculates Z-scores for depth data. It computes the median and Median Absolute Deviation (MAD) of
    the depth values across all contigs/regions and then calculates the Z-scores for each contig/region based on
    the MAD. Z-scores help identify peaks in the depth distribution.

    Example:
        depths_data = {
            'contig1': [10, 15, 12, 14],
            'contig2': [8, 7, 9, 10]
        }
        Z_scores, median, mad = get_ZScores(depths_data)
    """
    x = {}
    depths_list = np.empty(shape=0)

    # Combine depth values into a single numpy array
    for i in depths:
        depths_list = np.append(depths_list, depths[i])

    # Calculate median and MAD
    median = np.median(depths_list)
    mad = np.median(np.absolute(depths_list - median))




    # Check for MAD == 0 to avoid division by zero
    if mad == 0:
        if expect_mad_zero is True:
            plot_MAD_error_coverage(cov, output)
            logger.warning(
            "This error is likely because the wrong reads being mapped to genome because the MAD == 0.\nPlease check you are using correct reads.\nIf error persists please create an issue on GitHub (https://github.com/Chrisjrt/hafeZ)."
        )
            exit_error_gracefully_premap(output, start_time)
            sys.exit(0)
        else:
            exit_error_gracefully_premap(output, start_time)
            logger.error(
            "This error is likely because the wrong reads being mapped to genome because the MAD == 0.\n Use --expect_mad_zero to generate a MAD_ERROR_coverage_plot.\nPlease check you are using correct reads.\nIf error persists please create an issue on GitHub (https://github.com/Chrisjrt/hafeZ)."
        )

    # Calculate Z-scores for each contig/region
    # https://stats.stackexchange.com/questions/123895/mad-formula-for-outlier-detection
    for i in depths:
        x[i] = [((0.6745 * (x - median)) / mad) for x in depths[i]]
        x[i] = np.insert(
            x[i], 0, 0, axis=0
        )  # Add to start so peaks will find if the whole contig passes the threshold
        x[i] = np.append(x[i], 0)  # Add to the end for the same reason

    return x, median, mad


