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

out_pref = sys.argv[1]
output = f"{out_pref}/Quality/CheckM/checkm_cat.txt"

Path(output).touch()

# Open the output file for writing
with open(output, "w") as outfile:
    # Iterate through input files starting from the second argument
    for input_file in sys.argv[2:]:
        with open(input_file, "r") as infile:
            # Read and skip the header if not the first file
            if input_file != sys.argv[2]:
                next(infile)
            # Read and write the content of each file to the output
            for line in infile:
                outfile.write(line)