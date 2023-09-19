# ----------------------------------------------------------------------------
# Copyright (c) 2022--, Panpiper development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import getopt
import os
import numpy as np
import pandas as pd
import argparse
from collections import Counter

def checkm_filter(sample_list, file, completeness, contamination, N50):
    """
    Filter samples based off checkm parameters

    Parameters 
    ----------
    sample_list: string, required
    A text file which will contain names of all samples which pass the filters

    file: string, required
    A string of the filename of the checkm output, the log file

    completeness: int, default = 95
    An integer of the cutoff value at which to remove samples if they are under this number

    contamination: int, default = 5
    An integer of the cutoff value at which to remove samples if they are above this number

    N50: int, default = 5000
    An integer of the max N50 value, samples are filtered out if they exceed this value

    Returns
    ----------
    sample_pass: list
    The samples that passed the filters at this step returned as a list

    sample_fail: dataframe
    The samples that failed the filtes at this step returned 
    """

    # Reads the checkm text file
    df = pd.read_csv(file, sep='\t', header=0)
    # Includes completeness, contamination, and strain heterogeneity
    # Filter dataframe
    # For each row in the dataframe, only include a sample in a list
    # If it is above a certain value in each of the three categories
    # Completion above 95, Contamination below 5
    df_pass = df[df["Completeness"] >= completeness]
    df_pass = df_pass[df_pass["Contamination"] <= contamination]
    df_pass = df_pass[df_pass["Contig_N50"] >= N50]
    # Return the list of samples
    sample_pass = df_pass["Name"].tolist()
    sample_pass = [*set(sample_pass)]

    # Return the failed samples as a dataframe
    sample_fail = df[(df["Completeness"] < completeness) |
                     (df["Contamination"] > contamination) | 
                     (df_pass["Contig_N50"] < N50)]

    return sample_pass, sample_fail

def ani_filter(sample_list, file, ani, reference, output):
    """
    Filter samples based off checkm parameters

    Parameters 
    ----------
    sample_list: string, required
    A text file which will contain names of all samples which pass the filters

    file: string, required
    A string of the filename of the ani output, the matrix file

    ani: int, default = 95
    An integer of the cutoff value at which to remove samples if they are under this number compared to the reference

    reference: string, required
    The name of the reference sample as it appears in the ani matrix file

    Returns
    ----------
    sample_pass: list
    The samples that passed the filters at this step returned as a list

    sample_fail: dataframe
    The samples that failed the filtes at this step returned as a dataframe which contain information about why they failed
    """

    # Reads the ani text file
    df = pd.read_csv(file, sep='\t', header=0, index_col=None)
    df = df.fillna(value=0)
    df.columns = ["sample", "reference", "ani", "contig1", "contig2"]
    df["sample"] = df["sample"].apply(os.path.basename)
    df["reference"] = df["reference"].apply(os.path.basename)
    df["reference"] = df["reference"].str.replace(".fna", "")
    df["sample"] = df["sample"].str.replace(".fna", "")
    df.to_csv(output+'/fastani_reformat.csv', sep='\t')
    # Filter dataframe based off first column, select sample
    # Only add samples to the list if they are above 95% in the third column
    df_pass = df[df["ani"] >= ani]
    sample_fail = df.loc[df["ani"] < ani]
    # Add the samples from the second column to the list
    # Return the list of samples by APPENDING TO LIST
    sample_list = df_pass["sample"].tolist()
    sample_list = [*set(sample_list)]
    return sample_list, sample_fail


if __name__ == "__main__":
    checkm_list, ani_list = [], []
    # Result: 0 checkm,1 ani,2 outfile,3 completeness,4 contamination,
    #   5 strain_heterogeneity,6 ani_cutoff,7 reference
    # The default cut offs are based off the "Reasons an assembly is excluded from RefSeq" at the NCBI webserver
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", "--checkm")
    parser.add_argument("-a", "--ani")
    parser.add_argument("-o", "--outpath")
    parser.add_argument("-p", "--completeness", default=95, type=float)
    parser.add_argument("-c", "--contamination", default=5, type=float)
    parser.add_argument("-ac", "--ani_cuttoff", default=95, type=float)
    parser.add_argument("-r", "--reference")
    parser.add_argument("-n", "--N50", default=5000, type=float)
    args = parser.parse_args()
    print(args)
    if args.outpath:
        # create a file with the names of all the fastas in it
        if args.checkm:
            checkm_list, checkm_fail = checkm_filter(
                checkm_list, args.checkm, args.completeness, args.contamination, args.N50)
        if args.ani:
            ani_list, ani_fail = ani_filter(
                ani_list, args.ani, args.ani_cuttoff, args.reference, args.outpath)
        # Now I want to just keep the duplicates
        checkm_list += ani_list
        # Count the occurrences of each element in the combined list
        element_counts = Counter(checkm_list)
        # Filter the elements that occur exactly 2 times
        checkm_list = [element for element, count in element_counts.items() if count == 2]
        # merge dataframes of failed samples
        checkm_fail = checkm_fail.rename(columns={'Name': 'Sample'})
        if len(args.outpath) > 0:
            passed = args.outpath + "/sample_list.txt"
            with open(passed, mode='wt', encoding='utf-8') as myfile:
                myfile.write('\n'.join(checkm_list))
        checkm_fail.to_csv(args.outpath+'/failed_samples_checkm.csv', sep='\t')
        ani_fail.to_csv(args.outpath+'/failed_samples_ani.csv', sep='\t')