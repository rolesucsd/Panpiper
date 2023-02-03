'''============================================================================================
Renee Oles      10 June 2022
============================================================================================'''

import sys, getopt
import os
import numpy as np
import pandas as pd
import argparse


def checkm_filter(sample_list, file, completeness, contamination):
    """
    Filter samples based off checkm parameters
    """
    # Reads the checkm text file
    df = pd.read_csv(file, sep='\t', header=0)
    # Includes completeness, contamination, and strain heterogeneity
    # Filter dataframe
    # For each row in the dataframe, only include a sample in a list 
    # If it is above a certain value in each of the three categories
    # Completion above 90, Contamination below 20 or Contamination above 20 if Strain heterogenity below 20
    df_pass = df[df["Completeness"] >= completeness]
    df_pass = df_pass[df_pass["Contamination"] <= contamination] 
    # Return the list of samples
    sample_pass = df_pass["Bin_ID"].tolist()
    
    # Return the failed samples as a dataframe
    df_fail = df[(df["Completeness"] < completeness) | (df["Contamination"] > contamination)]
    return sample_pass, df_fail 
    #return sample_pass

if __name__ == "__main__":
    sample_list, checkm_list, ani_list = [], [], []
    # The default cut offs are based off the "Reasons an assembly is excluded from RefSeq" at the NCBI webserver
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--log")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-p", "--completeness", default=90, type=float)
    parser.add_argument("-c", "--contamination", default=20, type=float)
    args = parser.parse_args()
    print(args)
    if args.outfile:
        # create a file with the names of all the fastas in it
        if args.log:
            checkm_list, checkm_fail = checkm_filter(sample_list,args.log,args.completeness,args.contamination)
        # Now I want to just keep the duplicates
        checkm_fail = checkm_fail.rename(columns = {'Bin_ID':'Sample'})
        checkm_fail = pd.merge(checkm_fail, stats_fail, on = ['Sample'], how = 'outer')
        if len(args.outfile) > 0:
            with open(args.outfile, mode='wt', encoding='utf-8') as myfile:
                myfile.write('\n'.join(sample_list))
        checkm_fail.to_csv('metadata/to_bin.txt', sep ='\t')