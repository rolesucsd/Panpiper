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


def checkm_filter(sample_list, file, completeness, contamination):
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
    # Return the list of samples
    sample_pass = df_pass["Bin_ID"].tolist()
    sample_pass = [*set(sample_pass)]

    # Return the failed samples as a dataframe
    sample_fail = df[(df["Completeness"] < completeness) |
                     (df["Contamination"] > contamination)]
    return sample_pass, sample_fail


def checkm_stats_filter(sample_list, file, contig_number, N50):
    """
    Filter samples based off checkm parameters

    Parameters 
    ----------
    sample_list: string, required
    A text file which will contain names of all samples which pass the filters

    file: string, required
    A string of the filename of the checkm output, the stat file

    contig_number: int, default = 1000
    An integer of the cutoff value at which to remove samples if their assembly is greater than 1000 contigs

    N50: int, default = 5000
    An integer of the max N50 value, samples are filtered out if they exceed this value

    Returns
    ----------
    sample_pass: list
    The samples that passed the filters at this step returned as a list

    sample_fail: dataframe
    The samples that failed the filtes at this step returned as a dataframe which contain information about why they failed
    """

    # Reads the checkm text file
    df = pd.read_csv(file, sep='\t', header=0)
    # For each row in the dataframe, only include a sample in a list
    # If it is above a certain value in each of the categories
    # Contig num below input, N50 above 1e5
    # Subset the dataframe for the failed samples
    sample_fail = df[(df["# contigs"] > contig_number) | (df["N50 (contigs)"] < N50) ]
    df = df[df["# contigs"] <= contig_number]
    df = df[df["N50 (contigs)"] >= N50]
    # Return the list of samples and the the failed dataframe
    sample_list = df["Sample"].tolist()
    sample_list = [*set(sample_list)]
    return sample_list, sample_fail


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
    df.set_axis(["sample", "reference", "ani", "contig1", "contig2"],
                axis=1, inplace=True)
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
    sample_list, checkm_list, ani_list = [], [], []
    # Result: 0 checkm,1 ani,2 outfile,3 completeness,4 contamination,
    #   5 strain_heterogeneity,6 ani_cutoff,7 reference
    # The default cut offs are based off the "Reasons an assembly is excluded from RefSeq" at the NCBI webserver
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--log")
    parser.add_argument("-s", "--stat")
    parser.add_argument("-a", "--ani")
    parser.add_argument("-o", "--outpath")
    parser.add_argument("-p", "--completeness", default=95, type=float)
    parser.add_argument("-c", "--contamination", default=5, type=float)
    parser.add_argument("-sh", "--strain_heterogeneity", default=0, type=float)
    parser.add_argument("-ac", "--ani_cuttoff", default=95, type=float)
    parser.add_argument("-r", "--reference")
    parser.add_argument("-cn", "--contig_number", default=1000, type=float)
    parser.add_argument("-n", "--N50", default=5000, type=float)
    parser.add_argument("-lf", "--L50", default=500, type=float)
    args = parser.parse_args()
    print(args)
    if args.outpath:
        # create a file with the names of all the fastas in it
        if args.log:
            checkm_list, checkm_fail = checkm_filter(
                sample_list, args.log, args.completeness, args.contamination)
        if args.ani:
            ani_list, ani_fail = ani_filter(
                sample_list, args.ani, args.ani_cuttoff, args.reference, args.outpath)
        if args.stat:
            stats_list, stats_fail = checkm_stats_filter(
                sample_list, args.stat, args.contig_number, args.N50)
        # Now I want to just keep the duplicates
        checkm_list += ani_list
        checkm_list += stats_list
        uniques, uniques2 = [], []
        # remove all duplicates from list of samples that pass the filters
        for i in checkm_list:
            if i not in uniques:
                uniques += [i]
            elif i not in uniques2:
                uniques2 += [i]
            elif i not in sample_list:
                sample_list += [i]
        # merge dataframes of failed samples
        checkm_fail = checkm_fail.rename(columns={'Bin_ID': 'Sample'})
        checkm_fail = pd.merge(checkm_fail, stats_fail,
                               on=['Sample'], how='outer')
        if len(args.outpath) > 0:
            passed = args.outpath + "/sample_list.txt"
            with open(passed, mode='wt', encoding='utf-8') as myfile:
                myfile.write('\n'.join(sample_list))
        checkm_fail.to_csv(args.outpath+'/failed_samples_checkm.csv', sep='\t')
        ani_fail.to_csv(args.outpath+'/failed_samples_ani.csv', sep='\t')