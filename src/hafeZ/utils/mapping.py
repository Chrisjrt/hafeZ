from pathlib import Path
from hafeZ.utils.external_tools import ExternalTool

def minimap_paired(reads1: Path, reads2: Path, genome: Path, output: Path, threads: int, logdir: Path) -> None:
    """
    Align paired-end sequencing reads to a reference genome using Minimap2.

    Parameters:
        reads1 (Path): Path to the first set of input FASTQ files containing paired-end reads.
        reads2 (Path): Path to the second set of input FASTQ files containing paired-end reads.
        genome (Path): Path to the reference genome to align the reads against.
        output (Path): Path to the output directory where the alignment results will be stored.
        threads (int): Number of threads/cores to use for the alignment.
        logdir (Path): Path to the directory where log files will be stored.

    Returns:
        None

    This function uses Minimap2 to perform alignment of paired-end sequencing reads ('reads1' and 'reads2')
    against a reference genome ('genome'). The aligned output is stored in a SAM (Sequence Alignment/Map)
    format file named "temp_minimap.sam" in the 'output' directory.

    The 'threads' parameter specifies the number of CPU cores to use for the alignment, and 'logdir' is the
    directory where log files for the alignment process will be stored.
    """

    out_sam: Path = output / "temp_minimap.sam"

    minimap2 = ExternalTool(
        tool="minimap2",
        input=f"",
        output=f"",
        params=f'-ax sr {genome} {reads1} {reads2} -t {str(threads)} -o {out_sam} ',
        logdir=logdir,
    )
        
    ExternalTool.run_tool(minimap2)


def get_bam(output: Path, threads: int, memory_limit: str, logdir: Path) -> None:
    """
    Convert a SAM file to BAM format, then sort and index the resulting BAM file using Samtools.

    Parameters:
        output (Path): Path to the output directory where the BAM and sorted BAM files will be stored.
        threads (int): Number of threads/cores to use for the conversion and sorting.
        memory_limit (str): Memory allocation for sorting operations (e.g., '4G' for 4 gigabytes).
        logdir (Path): Path to the directory where log files will be stored.

    Returns:
        None

    This function takes a SAM format file and performs the following operations:
    
    1. Converts the SAM file to BAM format, storing it as "temp_minimap.bam" in the 'output' directory.
    2. Sorts the resulting BAM file using the specified number of threads and memory allocation,
       storing it as "temp_minimap_sorted.bam" in the 'output' directory.
    3. Indexes the sorted bam file

    The 'threads' parameter specifies the number of CPU cores to use for the conversion and sorting.
    The 'memory_limit' parameter sets the memory limit for sorting operations.
    """


    out_sam: Path = output / "temp_minimap.sam"
    out_bam: Path = output / "temp_minimap.bam"
    out_sorted_bam: Path = output / "temp_minimap_sorted.bam"

    samtools_view = ExternalTool(
        tool="samtools",
        input=f"",
        output=f"",
        params=f'view -S -b -@ {str(threads)} -o {out_bam} {out_sam} ',
        logdir=logdir,
    )

    ExternalTool.run_tool(samtools_view)

    samtools_sort = ExternalTool(
        tool="samtools",
        input=f"",
        output=f"",
        params=f'sort -@ {str(threads)} -m {memory_limit} -o {out_sorted_bam} {out_bam} ',
        logdir=logdir,
    )

    ExternalTool.run_tool(samtools_sort)

    samtools_index = ExternalTool(
        tool="samtools",
        input=f"",
        output=f"",
        params=f'index -@ {str(threads)} -b {out_sorted_bam}',
        logdir=logdir,
    )

    ExternalTool.run_tool(samtools_index)


def get_cov(output: Path, threads: int, logdir: Path) -> None: 
    """
    Calculate sequence coverage using Mosdepth on a sorted BAM file.

    Parameters:
        output (Path): Path to the output directory where coverage data will be stored.
        threads (int): Number of threads/cores to use for coverage calculation.
        logdir (Path): Path to the directory where log files will be stored.

    Returns:
        None

    This function calculates sequence coverage by using Mosdepth on a sorted BAM file. It generates coverage
    data for different regions of the genome and stores the results in the 'output' directory.

    The 'threads' parameter specifies the number of CPU cores to use for the coverage calculation.

    """

    mos_prefix: Path = output / "temp_minimap"
    out_sorted_bam: Path = output / "temp_minimap_sorted.bam"

    mosdepth = ExternalTool(
        tool="mosdepth",
        input=f"",
        output=f"",
        params=f'-t {str(threads)} -b 1 -n  {mos_prefix} {out_sorted_bam}',
        logdir=logdir,
    )

    ExternalTool.run_tool(mosdepth)

