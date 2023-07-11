import pandas as pd
import numpy as np
import re

def genes_long(genes, clusters, output):
    # Load data
    genes = pd.read_csv(genes, header=0, dtype={'Gene': str}, low_memory=False, comment='#')
    clusters = pd.read_csv(clusters, sep="\t")
    genes = genes.drop(genes.columns[[6, 7, 8, 9, 10]], axis=1)
    genes = genes.iloc[:, list([0, 2, 3]) + list(range(9, genes.shape[1]))]

    # Melt data
    genes_long = pd.melt(genes, id_vars=["Gene", "Annotation", "No. isolates"], var_name="Sample")
    genes_long = pd.merge(genes_long, clusters, on="Sample", how="left")
    genes_long = genes_long.rename(columns={genes_long.columns[4]: "contig"})

    genes_long['contig'] = genes_long['contig'].str.replace(r"[\[\]']", '', regex=True)
    # Split values in the 'contig' column by semicolon (;) and explode the resulting list into separate rows
    genes_long['contig'] = genes_long['contig'].str.split(';')
    genes_long = genes_long.explode('contig')

    # Save data to files
    genes_long.to_csv(output + "/genes_long.txt", sep="\t", index=False, quoting=None)

    # Create temporary dataframe
    genes_long['presence'] = 0
    genes_long.loc[genes_long['contig'].notna(), 'presence'] = 1

    return(genes_long) 

def gene_matrix_phylogroup(genes_long, output):
    # Group and summarize by Gene and Phylogroup
    genes_long = pd.read_csv(genes_long, sep="\t")
    genes_long['presence'] = 0
    genes_long.loc[genes_long['contig'].notna(), 'presence'] = 1

    # Group and summarize by Gene and Phylogroup
    genes_phylo = genes_long.iloc[:, [0, 3, 5, 6]]
    genes_phylo = genes_phylo.drop_duplicates()
    genes_phylo = genes_phylo.iloc[:, [0, 2, 3]]
    genes_phylo = genes_phylo.groupby(["Gene", "Phylogroup"]).agg({"presence": "sum"}).reset_index()


    # Pivot table and calculate percentages
    genes_phylo = genes_phylo.pivot(index='Gene', columns='Phylogroup', values='presence').fillna(0)
    genes_phylo.to_csv(output + "/genes_phylogroup_matrix.txt", sep="\t", index=True, quoting=None)

    genes_phylo_perc = round(genes_phylo.apply(lambda x: (x / max(x)) * 100, axis=0))
    # Rename columns and write to file
    genes_phylo_perc.reset_index(inplace=True)
    genes_phylo_perc.rename(columns={"Gene": "Gene"}, inplace=True)

    genes_phylo_perc.to_csv(output + "/genes_phylogroup_matrix_perc.txt", sep="\t", index=False, quoting=None)

def gene_matrix(genes,output):
    # Extract relevant columns
    genes = pd.read_csv(genes, header=0, dtype={'Gene': str}, low_memory=False, comment='#')
    genes = genes.drop(genes.columns[[6, 7, 8, 9, 10]], axis=1)
    genes_matrix = genes.iloc[:, [0, 9, *range(10, genes.shape[1])]]

    # Set gene name as index and remove gene column from matrix
    genes_matrix.index = genes_matrix['Gene']
    genes_matrix = genes_matrix.iloc[:, 1:]


    genes_matrix.to_csv(output + "/genes_matrix_names.txt", sep="\t", index=True, quoting=None)

    # Replace all NA values with 0
    genes_matrix = genes_matrix.fillna(0)    
    genes_matrix = genes_matrix.where(genes_matrix == 0, other=1)

    genes_matrix.to_csv(output + "/genes_matrix.txt", sep="\t", index=True, quoting=None)

    return(genes_matrix)

def process_data(bakta, genes, eggnog, output):
    # Load data
    bakta = pd.read_csv(bakta, delimiter="\t", header=0, comment='#')
    genes = pd.read_csv(genes, delimiter=",", low_memory=False, header=0, comment='#')
    eggnog = pd.read_csv(eggnog, delimiter="\t", header=0, comment='#')
    # Remove the '#' from the header row
    eggnog.columns = eggnog.columns.str.lstrip('#')

    # Clean data
    genes = genes.iloc[:, :6].join(genes.iloc[:, 11:], how='outer')
    bakta = bakta.rename(columns={bakta.columns[2]:"Gene_Name" , bakta.columns[0] : "Gene"})
    bakta = genes.iloc[:, :8].merge(bakta, on="Gene", how="left")
    bakta.drop_duplicates(subset="Gene", keep="first", inplace=True)
    eggnog = eggnog.rename(columns={eggnog.columns[0]:"Gene"})
    genes_anno = bakta.merge(eggnog, on="Gene", how="left").drop_duplicates(subset="Gene", keep="first")
    genes_anno["group"] = "ACCESSORY"
    genes_anno.loc[genes_anno["No. isolates"] >= (genes_anno["No. isolates"].max() * 0.99), "group"] = "CORE"
    genes_anno.loc[genes_anno["No. isolates"] / genes_anno["No. isolates"].max() <= 0.01, "group"] = "UNIQUE"
    gene_columns = list(genes.columns)

    #Write output
    genes_anno.to_csv(f"{output}/genes_anno.txt", sep='\t', index=False)

    return(genes_anno)