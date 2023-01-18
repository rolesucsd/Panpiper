'''============================================================================================
Renee Oles      4 Oct 2021
============================================================================================'''

import sys
import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage

def file_to_heatmap(file, output, path=""):
    """
    Open the matrix ani text file and make a heatmap
    """
    # Reads the ani text file
    df = pd.read_table(file, sep='\t', index_col=0, header=0)
    # Make heatmap
    df.columns = df.columns.str.strip(path)
    df.index = df.index.str.strip(path)
    # Fill empty cells with 0 to avoid errors
    df = df.fillna(0)
    # calculate minimum for plot
    array = df.to_numpy()
    minim = array.min()
    link = linkage(df) # D being the measurement
    # color map
    plt.figure(dpi=1200)
    sb.clustermap(df, annot=False, vmin=minim, vmax=100, 
        cmap="PuOr_r", square=True, row_linkage=link, col_linkage=link)
    print(sb.color_palette("PuOr_r").as_hex())
    # title
    title = 'Average Nucleotide Identity'
    plt.title(title, loc='left', fontsize=12)
    plt.axis('off')
    plt.autoscale()
    plt.show()
    plt.savefig(output, bbox_inches="tight")
    return 0

if __name__ == "__main__":
    if len(sys.argv) == 3:
        file_to_heatmap(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        file_to_heatmap(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        print("Enter file, output name, and optional path to remove from header")