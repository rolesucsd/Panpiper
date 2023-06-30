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

args = sys.argv[1:]
out_pref = args[-1]
output = f"{out_pref}/Quality/CheckM/checkm_log.txt"
args = args[:-1]

Path(output).touch()

def process_file(input_file):
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
                            with open(output, 'a') as out_file:
                                out_file.write(line + "\n")
    except FileNotFoundError:
        # If the file doesn't exist, skip the normal part and fill the line with NAs
        line = "\t".join(["NA"] * 5) + "\n"
        with open(output, 'a') as out_file:
            out_file.write(line)

for file_path in args:
    process_file(file_path)

# LOAD DATA
df_c = pd.DataFrame()
with open(output, 'r') as file:
    df_c = pd.read_table(file, header=None, sep="\t")

df_c.rename(columns={0: "first"}, inplace=True)
df_c = df_c["first"].str.split("\s{2,}", expand=True)
print(df_c)
df_c = df_c[[0, 1, 11, 12, 13]]
df_c.columns = ["Bin_ID", "Marker_lineage", "Completeness", "Contamination", "Strain_heterogeneity"]
#df_c[["Blank", "Bin_ID", "Marker_lineage", "Marker_lineage_ID", "genomes", "markers", "marker_sets", "0", "1", "2", "3", "4", "5+", "Completeness", "Contamination", "Strain_heterogeneity"]] = 
print(df_c)
df_c["Contamination"] = pd.to_numeric(df_c["Contamination"])
df_c["Completeness"] = pd.to_numeric(df_c["Completeness"])
df_c["Strain_heterogeneity"] = pd.to_numeric(df_c["Strain_heterogeneity"])
df_c.loc[(df_c["Strain_heterogeneity"] == 100) & (df_c["Contamination"] == 100), "Strain_heterogeneity"] = 0
df_c.to_csv(output, sep="\t", index=False)

filtered_df = df_c[(df_c["Contamination"] < 5) & (df_c["Completeness"] > 95)]
filtered_df.iloc[:, 0].to_csv(f"{out_pref}/Quality/CheckM/checkm_log_filter.txt", sep="\t", index=False, header=False)
