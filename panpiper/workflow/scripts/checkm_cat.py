# ----------------------------------------------------------------------------
# Copyright (c) 2022--, Panpiper development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd
import re
import sys
from pathlib import Path

def concatenate_checkm_files(out_pref: str, input_files: list):
    """
    Concatenate content of multiple CheckM output files into a single file.

    Args:
        out_pref (str): Output prefix.
        input_files (list): List of input file paths.

    Returns:
        None
    """
    output = f"{out_pref}/Quality/CheckM/checkm_cat.txt"

    Path(output).touch()

    # Open the output file for writing
    with open(output, "w") as outfile:
        # Iterate through input files starting from the second argument
        for input_file in input_files:
            with open(input_file, "r") as infile:
                # Read and skip the header if not the first file
                if input_file != input_files[0]:
                    next(infile)
                # Read and write the content of each file to the output
                for line in infile:
                    outfile.write(line)

if __name__ == "__main__":
    # Get output prefix and input file paths from command-line arguments
    out_pref = sys.argv[1]
    input_files = sys.argv[2:]

    # Call the function
    concatenate_checkm_files(out_pref, input_files)
