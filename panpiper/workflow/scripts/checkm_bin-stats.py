import sys
import os
import pandas as pd

def process_checkm_files(output_prefix: str, input_files: list):
    """
    Process CheckM output files and generate a consolidated statistics file.

    Args:
        output_prefix (str): Prefix for the output file path.
        input_files (list): List of input file paths.

    Returns:
        None
    """
    # Define output path
    output_path = os.path.join(output_prefix, "Quality/CheckM/checkm_stats.txt")

    checkm = pd.DataFrame()

    # Read and process input files
    for input_file in input_files:
        line = pd.read_csv(input_file, sep='\t', header=None)
        checkm = pd.concat([checkm, line], ignore_index=True)

    # Process the checkm2 column
    checkm2 = checkm.drop(checkm.columns[0], axis=1)
    checkm2 = checkm2.astype(str).applymap(lambda x: str(x).replace("'", "").replace("{", "").replace("}", ""))
    checkm2 = checkm2[1].str.split(',', expand=True)
    checkm2 = checkm2.applymap(lambda x: x.split(":")[1])
    checkm2.columns = ["GC", "GC std", "Genome size", "ambiguous bases", "scaffolds", "contigs",
                       "Longest scaffold", "Longest contig", "N50 scaffolds", "N50 contigs", "Mean scaffold length",
                       "Mean contig length", "Coding density", "Translation table", "predicted genes"]

    # Process the Sample column
    checkm2.insert(0, "Sample", checkm[0])

    # Save the DataFrame to a tab-separated text file
    checkm2.to_csv(output_path, sep='\t', index=False)

# Get input file paths from command-line arguments
input_files = sys.argv[2:]

# Define output path
output_prefix = sys.argv[1]

# Call the function
process_checkm_files(output_prefix, input_files)
