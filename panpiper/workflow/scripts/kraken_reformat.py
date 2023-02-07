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


def agg_df(x):
    return pd.Series(dict(Length=x['Length'].sum(),
                          Contig="%s" % ', '.join(x['Contig'])))


def reformat_kraken(basename, krakenpath):
    """
    Changes the kraken output to be a dataframe 
    Parameters 
    ----------
    basename: string, required
        The name of the sample and basename of the file

    krakenpath: string, required
        The name of out directory, and the directory where the kraken files are stored
    """
    # Read in the file
    filename = krakenpath+"/"+basename+".out"

    df = pd.read_csv(filename, sep='\t', header=None, names=[
                     "D1", "Contig", "ID", "Length", "D2"])

    # Remove unnecessary columns
    df = df[["ID", "Length", "Contig"]]

    # Get the sum of the column, Length
    df_sum = df['Length'].sum()

    # Aggregate by the second column, concatenate all the contig IDs together in a list
    df = df.groupby('ID').agg({'Length': 'sum', 'Contig': ','.join})

    # add a column called percent which is the percent of the sum for each line
    df["Percent"] = df.Length / df_sum * 100

    # add another column which is the base name this column should be called
    df["Sample"] = basename

    # reformat again to : Sample, ID, Percent, Length, Contig
    df = df[['Sample', 'Percent', 'Length', 'Contig']]

    df.to_csv(krakenpath+"/"+basename+"_reformat.txt", header=True, sep='\t')


if __name__ == "__main__":
    basename = sys.argv[1]
    krakenpath = sys.argv[2]
    reformat_kraken(basename, krakenpath)
