#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --ntasks=5
#SBATCH --mem=100G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roles@health.ucsd.edu

cd /panfs/roles/Panpiper

source /home/roles/anaconda3/bin/activate

#conda create -n panpiper -c bioconda -c conda-forge snakemake mamba 

conda activate panpiper

#bakta_db download --output /panfs/roles/Panpiper/panpiper/databases/bakta

#panpiper -o /panfs/roles/BF --unlock 
#panpiper -o /panfs/roles/BF -q /panfs/roles/Isolates/Reference/All --log_lvl DEBUG -w assembly --cluster_type slurm --cluster_config cluster.json --cluster_args "sbatch -A {cluster.account} --mem {cluster.mem} -t {cluster.time} --cpus-per-task {cluster.cpus}"

panpiper -o /panfs/roles/BF --log_lvl DEBUG -a /panfs/roles/BF/Assembly -s  /panfs/roles/BF/Assembly/complete_assembly_files.txt -r /panfs/roles/BF/reference/9343.fna -w quality --cluster_type slurm --cluster_config cluster.json --cluster_args "sbatch -A {cluster.account} --mem {cluster.mem} -t {cluster.time} --cpus-per-task {cluster.cpus}"
#panpiper -o /panfs/roles/BF --log_lvl DEBUG -a /panfs/roles/BF/Quality/Assembly_filter -s  /panfs/roles/BF/Quality/sample_list.txt -r /panfs/roles/Package_test/reference/9343.fna -w pangenome --cluster_type slurm --cluster_config cluster.json --cluster_args "sbatch -A {cluster.account} --mem {cluster.mem} -t {cluster.time} --cpus-per-task {cluster.cpus}"

# TODO: Assembly - some assemblies get stuck in shovill
# TODO: Assembly - can't find move_files.sh???? 
# TODO: Quality - I don't include completeness and contamination at the moment 
