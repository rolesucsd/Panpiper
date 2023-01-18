#============================================================================================
#Renee Oles      9 Nov 2021
#============================================================================================

import sys
import os
import numpy as np
import pandas as pd

def amr_concat(files, output, path=""):
    """
    Concatenate the output of amr 
    """
    output = path+output
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

if __name__ == "__main__":
    print(sys.argv)
    amr_concat(sys.argv, "AMR.txt", "../MGWAS/")