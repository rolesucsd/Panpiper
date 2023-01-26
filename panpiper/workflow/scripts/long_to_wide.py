# ----------------------------------------------------------------------------
# Copyright (c) 2022--, Panpiper development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import os
import numpy as np
import pandas as pd

def long_to_wide(file, output, path=""):
    """
    Open the text file and convert format from long to wide (matrix). Writes the output to the file specified.

    Parameters 
    ----------
    file: string, required
        The filename of the file read in 

    output: string, required
        The filename of the desired output

    path: string
        The full path of the file

    """
    # Reads the ani text file
    df = pd.read_csv(file, sep='\t', names=["Reference","Sample","PercentIdentity","Contig1","Contig2"])
    # TODO: Change names - remove path 
    for row in range(0,len(df.index)):
        df.at[row,'Reference'] = os.path.basename(df.at[row,'Reference'])
        df.at[row,'Sample'] = os.path.basename(df.at[row,'Sample'])
    df['Reference'] = df['Reference'].str.replace('\\.fasta', '')
    df['Reference'] = df['Reference'].str.replace('\\.fna', '')
    df['Sample'] = df['Sample'].str.replace('\\.fasta', '')
    df['Sample'] = df['Sample'].str.replace('\\.fna', '')
    df = df.pivot_table(index='Reference', columns='Sample', values='PercentIdentity', aggfunc='mean')
    if os.path.exists(output):
        os.remove(output)
        print("File exists, overwriting")
    df.to_csv(output, header=True, index=True, sep='\t', mode='a')

    return 0

if __name__ == "__main__":
    if len(sys.argv) == 3:
        long_to_wide(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        long_to_wide(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        print("Enter file, output name, and optional path to remove from header")