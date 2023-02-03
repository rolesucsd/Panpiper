#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --ntasks=5
#SBATCH --mem=200G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roles@health.ucsd.edu

cd /panfs/roles/Panpiper

source /home/roles/anaconda3/bin/activate

#conda create -n panpiper -c bioconda -c conda-forge snakemake mamba 

conda activate fasttree

FastTree -gtr -nt -gamma -seed 12345 /panfs/roles/BF/Pangenome/Panaroo/core_gene_alignment.aln > /panfs/roles/BF/Pangenome/Phylogeny/fasttree.nwk
#bakta_proteins --db panpiper/databases/bakta --output /panfs/roles/BF/Pangenome/Summary/ --prefix pan_genome_reference /panfs/roles/BF/Pangenome/Panaroo/pan_genome_reference.faa &> /panfs/roles/BF/report/bakta_pan.log