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

def touch_file(file_path):
    """
    Create an empty file if it doesn't exist.

    Args:
        file_path (str): Path to the file.

    Returns:
        None
    """
    Path(file_path).touch()

def process_file(input_file, output_file):
    """
    Process a CheckM log file and write the results to an output file.

    Args:
        input_file (str): Path to the input file.
        output_file (str): Path to the output file.

    Returns:
        None
    """
    try:
        with open(input_file, 'r') as file:
            for line in file:
                if re.search("---", line):
                    while True:
                        line = file.readline().strip()
                        if not line:
                            break
                        if not re.search("Bin", line) and not re.search("---", line) and not re.search("\\[", line):
                            input_file = input_file.replace("*Quality/Assembly_filter/", "")
                            input_file = input_file.replace("/lineage.log", "")
                            line = line.replace("filter", input_file)
                            with open(output_file, 'a') as out_file:
                                out_file.write(line + "\n")
    except FileNotFoundError:
        # If the file doesn't exist, skip the normal part and fill the line with NAs
        line = "\t".join(["NA"] * 5) + "\n"
        with open(output_file, 'a') as out_file:
            out_file.write(line)

def load_data_and_process(output_file):
    """
    Load data from a processed file and perform additional processing.

    Args:
        output_file (str): Path to the output file.

    Returns:
        pd.DataFrame: Processed DataFrame.
    """
    df_c = pd.DataFrame()
    with open(output_file, 'r') as file:
        df_c = pd.read_table(file, header=None, sep="\t")

    df_c.rename(columns={0: "first"}, inplace=True)
    df_c = df_c["first"].str.split("\s{2,}", expand=True)
    df_c = df_c[[0, 1, 11, 12, 13]]
    df_c.columns = ["Bin_ID", "Marker_lineage", "Completeness", "Contamination", "Strain_heterogeneity"]

    df_c["Contamination"] = pd.to_numeric(df_c["Contamination"])
    df_c["Completeness"] = pd.to_numeric(df_c["Completeness"])
    df_c["Strain_heterogeneity"] = pd.to_numeric(df_c["Strain_heterogeneity"])
    df_c.loc[(df_c["Strain_heterogeneity"] == 100) & (df_c["Contamination"] == 100), "Strain_heterogeneity"] = 0

    return df_c

def filter_and_save_results(df, out_pref):
    """
    Filter DataFrame and save results to a file.

    Args:
        df (pd.DataFrame): Input DataFrame.
        out_pref (str): Output prefix.

    Returns:
        None
    """
    filtered_df = df[(df["Contamination"] < 5) & (df["Completeness"] > 95)]
    filtered_df.iloc[:, 0].to_csv(f"{out_pref}/Quality/CheckM/checkm_log_filter.txt", sep="\t", index=False, header=False)

def main():
    args = sys.argv[1:]
    out_pref = args[-1]
    output_file = f"{out_pref}/Quality/CheckM/checkm_log.txt"
    args = args[:-1]

    touch_file(output_file)

    for file_path in args:
        process_file(file_path, output_file)

    processed_df = load_data_and_process(output_file)
    filter_and_save_results(processed_df, out_pref)

if __name__ == "__main__":
    main()
