'''============================================================================================
Renee Oles      12 Nov 2021
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
    sample_pass = [*set(sample_pass)]
    
    # Return the failed samples as a dataframe
    df_fail = df[(df["Completeness"] < completeness) | (df["Contamination"] > contamination)]
    return sample_pass, df_fail 
    #return sample_pass

def checkm_stats_filter(sample_list, file, genome_size, contig_number, N50, GC):
    """
    Filter samples based off checkm bin stats parameters
    """
    # Reads the checkm text file
    df = pd.read_csv(file, sep='\t', header=0)
    # For each row in the dataframe, only include a sample in a list 
    # If it is above a certain value in each of the categories
    # Genome size within 50% of ref, Contig num below input, N50 above 1e5, GC within 30% of ref
    mini = genome_size-genome_size*.5
    maxi = genome_size+genome_size*.5
    print(mini, maxi)
    # Subset the dataframe for the failed samples
    df_fail = df[(df["contigs"] > contig_number) | (df["N50 contigs"] < N50) | (df["Genome size"] < int(mini)) | 
        (df["Genome size"] > int(maxi)) | (df["GC"] > (GC+GC*.1)/100) | (df["GC"] < (GC-GC*.1)/100)]

    df = df[df["contigs"] <= contig_number]
    df = df[df["N50 contigs"] >= N50]
    if genome_size != 0:
        df = df[df["Genome size"] >= int(mini)] 
        df = df[df["Genome size"] <= int(maxi)]
    df = df[df["GC"] <= (GC+GC*.1)/100]
    df = df[df["GC"] >= (GC-GC*.1)/100]
    # Return the list of samples and the the failed dataframe
    sample_list = df["Sample"].tolist()
    sample_list = [*set(sample_list)]
    return sample_list, df_fail

def ani_filter(sample_list, file, ani, reference):
    """
    Filter samples based off ani parameters
    """
    # Reads the ani text file
    df = pd.read_csv(file, sep='\t', header=0, index_col=0)
    df = df.fillna(value = 0)
    # Filter dataframe based off first column, select sample
    # Only add samples to the list if they are above 95% in the third column 
    df_pass = df[df[reference] >= ani]
    df_fail = df.loc[df[reference] < ani, reference]
    # Add the samples from the second column to the list
    # Return the list of samples by APPENDING TO LIST 
    sample_list = list(df_pass.index)
    sample_list = [*set(sample_list)]
    return sample_list, df_fail

if __name__ == "__main__":
    sample_list, checkm_list, ani_list = [], [], []
    # Result: 0 checkm,1 ani,2 outfile,3 completeness,4 contamination,
    #   5 strain_heterogeneity,6 ani_cutoff,7 reference
    # The default cut offs are based off the "Reasons an assembly is excluded from RefSeq" at the NCBI webserver
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--log")
    parser.add_argument("-s", "--stat")
    parser.add_argument("-a", "--ani")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-p", "--completeness", default=95, type=float)
    parser.add_argument("-c", "--contamination", default=5, type=float)
    parser.add_argument("-sh", "--strain_heterogeneity", default=0, type=float)
    parser.add_argument("-f", "--ani_cuttoff", default=85, type=float)
    parser.add_argument("-r", "--reference")
    parser.add_argument("-g", "--genome_size", default=0, type=float)
    parser.add_argument("-cn", "--contig_number", default=1000, type=float)
    parser.add_argument("-n", "--N50", default=5000, type=float)
    parser.add_argument("-lf", "--L50", default=500, type=float)
    parser.add_argument("-gc", "--GC", default=0, type=float)
    args = parser.parse_args()
    print(args)
    if args.outfile:
        # create a file with the names of all the fastas in it
        if args.log:
            checkm_list, checkm_fail = checkm_filter(sample_list,args.log,args.completeness,args.contamination)
        if args.ani:
            ani_list, ani_fail = ani_filter(sample_list,args.ani,args.ani_cuttoff,args.reference)
        if args.stat:
            stats_list, stats_fail = checkm_stats_filter(sample_list,args.stat,args.genome_size,args.contig_number,args.N50,args.GC)
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
        checkm_fail = checkm_fail.rename(columns = {'Bin_ID':'Sample'})
        checkm_fail = pd.merge(checkm_fail, stats_fail, on = ['Sample'], how = 'outer')
        if len(args.outfile) > 0:
            with open(args.outfile, mode='wt', encoding='utf-8') as myfile:
                myfile.write('\n'.join(sample_list))
        checkm_fail.to_csv('logs/failed_samples_checkm.csv', sep ='\t')
        ani_fail.to_csv('logs/failed_samples_ani.csv', sep ='\t')
        