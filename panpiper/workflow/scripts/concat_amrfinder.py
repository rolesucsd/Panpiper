#============================================================================================
#Renee Oles      9 Nov 2021
#============================================================================================

import sys, getopt
import os
import numpy as np
import pandas as pd
import argparse

def amr_concat(files, output):
    """
    Concatenate the output of amr 
    """

    df = pd.DataFrame()
    for f in files:
        df_new = pd.read_csv(f, sep="\t")
        df_new['Sample'] = f
        df = pd.concat([df_new,df])

    if os.path.exists(output):
        os.remove(output)
        print("File exists, overwriting")
    df.to_csv(output, header=True, index=False, sep='\t', mode='a')
    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', nargs='+', help='List of input files')
    parser.add_argument('-o', '--output', help='Output file name')
    args = parser.parse_args()

    amr_concat(args.input, args.output)