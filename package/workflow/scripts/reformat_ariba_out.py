'''============================================================================================
Renee Oles      9 Nov 2021
============================================================================================'''

import sys
import os
import numpy as np
import pandas as pd

def csv_to_tab(file, output, path=""):
    """
    Open the text file and convert format)
    """
    # Reads the ani text file
    df = pd.read_csv(file)
    #Change names of samples by removing the path
    df['name'] = df['name'].str.replace('/report.tsv', '')
    df['name'] = df['name'].str.replace('../MGWAS/Ariba/', '')
    
    # Change yes/no of columns to 1 0 
    for column in df:
        if column != 'name':
            df[column] = df[column].map({'yes': 1, 'no': 0})

    if os.path.exists(output):
        os.remove(output)
        print("File exists, overwriting")
    df.to_csv(output, header=True, index=False, sep='\t', mode='a')
    return 0

if __name__ == "__main__":
    if len(sys.argv) == 3:
        csv_to_tab(sys.argv[1], sys.argv[2])
    else:
        print("Enter file and output name")