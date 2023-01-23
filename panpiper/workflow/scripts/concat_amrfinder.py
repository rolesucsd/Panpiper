import sys, getopt
import os
import numpy as np
import pandas as pd
import argparse

def amr_concat(files, output):
    """
    Concatenate the output of AMR into master AMR matrix 

    Parameters 
    ----------
    files: list of strings, required
    A list of input files from the output of AMRFinderPlus

    output: string, required
    A string representing the desired output filename, should include path

    Returns
    ----------
    Writes the matrix to the output file. 

    TODO: Should raise exceptions or write warnings
    """

    df = pd.DataFrame()
    # For each filename in the list, read it in and concat it to the new dataframe, df
    for f in files:
        df_new = pd.read_csv(f, sep="\t")
        df_new['Sample'] = f
        df = pd.concat([df_new,df])

    # If the output file exists, overwrite it
    if os.path.exists(output):
        os.remove(output)
        print("File exists, overwriting")

    # Write the dataframe to the desired location
    df.to_csv(output, header=True, index=False, sep='\t', mode='a')
    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', nargs='+', help='List of input files')
    parser.add_argument('-o', '--output', help='Output file name')
    args = parser.parse_args()

    amr_concat(args.input, args.output)