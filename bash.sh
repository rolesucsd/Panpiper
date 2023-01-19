#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --ntasks=5
#SBATCH --mem=100G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roles@health.ucsd.edu

cd /panfs/roles/Panpiper/panpiper

source /home/roles/anaconda3/bin/activate

#conda create -n panpiper -c bioconda -c conda-forge snakemake mamba 

conda activate panpiper

python main.py -o ../../Package_test --unlock 
#python main.py -o ../../Package_test -q ../../Package_test/fastq -w assembly
#python main.py -o ../../Package_test -a ../../Package_test/Assembly -s  ../../Package_test/Assembly/complete_assembly_files.txt -r /panfs/roles/Package_test/reference/9343.fna -w quality

python main.py -o ../../Package_test -a ../../Package_test/Quality/Assembly_filter -s  ../../Package_test/Quality/sample_list.txt -r /panfs/roles/Package_test/reference/9343.fna -w pangenome --cluster_type slurm --cluster_config cluster.json --cluster_args "sbatch -A {cluster.account} --mem {cluster.mem} -t {cluster.time} --cpus-per-task {cluster.cpus}"