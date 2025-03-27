import pandas as pd
import numpy as np
import pathlib
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import pysam
import os
import argparse
import pyfaidx
import sys
import re

##### helper functions for processing input BAM files
def fetch_reads(bamfile, region):
    """
    Fetches reads from a BAM file at a specific region.
    Parameters:
    - bamfile (str): Path to the BAM file.
    - region (str): Genomic region to fetch reads from. This should be taken from the file CpG_clusters_whole_genome_radius_100.bed
    Returns:
    - readdf (pandas.DataFrame): DataFrame containing the fetched reads with the following columns:
        - chrom: Chromosome name
        - start: Start position of the read
        - cigar: CIGAR string
        - flen: Length of the read
        - seq: Sequence of the read
        - methyl_string: Methyl string
        - XR: XR tag value
        - XG: XG tag value
        - sample: Sample name
        - region: Genomic region
    Example: 
    bamfile = "/Users/hieunguyen/data/bam_files/highdepth_cancer_WGBS_bismark.bam"
    region = "1:1103264-1103363"
    """
    # Function implementation goes here
    pass
    
    bamfile_obj = pysam.AlignmentFile(bamfile).fetch(region = region)
    reads = []
    for read in bamfile_obj:
        reads.append(read)
    readdf = pd.DataFrame()
    readdf["chrom"] = [read.to_dict()["ref_name"] for read in reads]
    readdf["start"] = [read.to_dict()["ref_pos"] for read in reads]
    readdf["cigar"] = [read.to_dict()["cigar"] for read in reads]
    readdf["flen"] = [read.to_dict()["length"] for read in reads]
    readdf["seq"] = [read.to_dict()["seq"] for read in reads]
    # readdf["methyl_string"] = [read.to_dict()["tags"][2] for read in reads]
    # readdf["XR"] = [read.to_dict()["tags"][3] for read in reads]
    # readdf["XG"] = [read.to_dict()["tags"][4] for read in reads]
    readdf["sample"] = str(bamfile).split("/")[-1].split(".")[0]
    readdf["region"] = region
    return readdf

##### helper functions for processing CpG status in reads#
def get_refseq(path_to_all_fa, chrom, start, end):
    """
    Retrieves the reference sequence from a given FASTA file.
    Args:
        path_to_all_fa (str): The path to the directory containing all the FASTA files.
        chrom (str): The chromosome identifier.
        start (int): The starting position of the sequence.
        end (int): The ending position of the sequence.
    Returns:
        str: The uppercase reference sequence.
    Raises:
        FileNotFoundError: If the FASTA file for the specified chromosome is not found.
    """
    refseq = pyfaidx.Fasta(os.path.join(path_to_all_fa, "{}.fa".format(chrom)))
    return(str.upper(refseq.get_seq(name = "{}".format(chrom), start = start, end = end).seq))
            


###### calculate the distance between the read and a nearest CpG site, if exists
def get_min_dist_to_cpg(chrom, start, path_to_all_fa, radius = 5):
    refseq_at_cluster = get_refseq(path_to_all_fa = path_to_all_fa, 
                                    chrom = chrom, 
                                    start = start - radius - 1, 
                                    end = start + radius)
    cpg_dists = [m.start(0) - radius for m in re.finditer("CG", refseq_at_cluster)]
    if len(cpg_dists) == 0:
        min_dist_to_cpg = "NA"
    else:
        min_abs_dist = min([abs(item) for item in cpg_dists])
        min_dist_to_cpg = [item  for item in cpg_dists if abs(item) == min_abs_dist][0]
    return min_dist_to_cpg
