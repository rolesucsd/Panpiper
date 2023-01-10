#!/usr/bin/env python

import argparse
import logging
import os
import sys
import pkg_resources

from os.path import join, dirname, basename
from shutil import copyfile
from argparse import RawTextHelpFormatter

from package.controller import Controller
from package.workflow import Workflow

# Snakefiles
_ROOT = os.path.abspath(os.path.dirname(__file__))
WORKFLOW_ASSEMBLY = os.path.join(_ROOT, 'workflow', 'assembly.smk')
WORKFLOW_QUALITY = os.path.join(_ROOT, 'workflow', 'quality.smk')
WORKFLOW_PANGENOME = os.path.join(_ROOT, 'workflow', 'pangenome.smk')

def cli():
    
    ########## Arguments ##########
    ap = argparse.ArgumentParser(description='Package version {}'.format(pkg_resources.require("package")[0].version), add_help=False, allow_abbrev=False, formatter_class=RawTextHelpFormatter)
    
    # Required
    apr = ap.add_argument_group('main arguments')
    apr.add_argument('-q', '--fastq', help='Path to fastq file directory', required=False, default="skip", type=str)
    apr.add_argument('-a', '--fasta', help='Path to fasta file directory', required=False, default="skip", type=str)
    apr.add_argument('-s', '--sample_list', help='Line delimited file of samples passing quality control', required=False, default="skip", type=str)
    apr.add_argument('-o', '--output', help='Prefix for output directory', required=True)
    apr.add_argument('-r', '--reference', help='Reference fasta file, includes path', type=str, required=False, default="skip", type=str)

    # Cluster arguments
    apc = ap.add_argument_group('compute cluster arguments')
    apc.add_argument('--cluster', help='Cluster compute structure [%(default)s]', default=None, type=str, choices=[None,'qsub','slurm','drmaa'])
    apc.add_argument('--cluster_info', help='Cluster scheduler arguments when submitting cluster jobs.\nHas to contain the following special strings:\n{memory}, {cores}, and {runtime}.\nThese special strings will be substituted by maginator to indicate resources for each job.\n{memory} is substituted for the memory in GB.\n{runtime} is substituted with the time in the following format: DD:HH:MM:SS.\nCan also contain user names, groups, etc. required by the cluster scheduler', default=None, type=str)
    apc.add_argument('--max_jobs', help='Maximum number of cluster jobs [%(default)s]', default=50, type=int)
    
    # Optional
    apo = ap.add_argument_group('optional arguments')
    apo.add_argument("-h", "--help", action="help", help="show this help message and exit")
    apo.add_argument('--max_cores', help='Maximum number of cores [%(default)s]', default=40, type=int)
    apo.add_argument('--max_mem', help='Maximum memory in GB [%(default)s]', default=180, type=int)
    apo.add_argument('--log_lvl', help='Logging level [%(default)s].', default='INFO', type=str, choices=['DEBUG','INFO','WARNING','ERROR'])
    apo.add_argument('--only_conda', help='Only install conda environments, then exit', action='store_true')
    apo.add_argument('--snake', help='Only run specific snakemake command. For debug purposes', type=str)
    apo.add_argument('--unlock', help='Unlock snakemake directory in case of unexpected exists, then exit', action='store_true')

    # Parameters
    # TODO Change parameters
    app = ap.add_argument_group('parameters')
    app.add_argument('--gc', help='Minimum bin size for inclusion [%(default)s].', default=200000, type=int)
    app.add_argument('--genome_size', help='If set, bin clusters will not be aggregated to metagenomic species, but treated separately.', action='store_true')
    app.add_argument('--annotation_prevalence', help='Minimum prevalence of taxonomic assignment in a cluster of bins to call consensus [%(default)s]', default=0.75, type=float)
    app.add_argument('--clustering_coverage', help='Alignment coverage for clustering of genes with MMseqs2 [%(default)s]', default=0.95, type=float)
    app.add_argument('--clustering_min_seq_id', help='Sequence identity threshold for clustering of genes with MMseqs2 [%(default)s]', default=0.95, type=float)
    app.add_argument('--clustering_type', help='Sequence type for gene clustering with MMseqs2. Nucleotide- or protein-level [%(default)s]', default='protein', type=str, choices=['nucleotide', 'protein'])
    app.add_argument('--min_gtdb_markers', help='Minimum GTDBtk marker genes shared between MGS and outgroup for rooting trees [%(default)s]', default=10, type=int)
    app.add_argument('--marker_gene_cluster_prevalence', help='Minimum prevalence of marker genes to be selected for rooting of MGS trees [%(default)s]', default=0.5, type=float)
    app.add_argument('--min_af', help='Minimim allele frequency for calling a base when creating phylogenies [%(default)s]', default=0.8, type=float)
    app.add_argument('--min_depth', help='Minimim read depth for calling a base when creating phylogenies [%(default)s]', default=2, type=int)
    app.add_argument('--min_nonN', help='Minimum fraction of non-N characters of a sample to be included in a phylogeny [%(default)s]', default=0.5, type=float)
    app.add_argument('--min_marker_genes', help='Minimum marker genes to be detected for inclusion of a sample in a phylogeny [%(default)s]', default=2, type=int)
    app.add_argument('--min_signature_genes', help='Minimum signature genes to be detected for inclusion of a sample in a phylogeny [%(default)s]', default=50, type=int)
    app.add_argument('--phylo', help='Software for phylogeny inference. Either fast (fasttree) or slow and more accurate (iqtree) [%(default)s]', default='fasttree', type=str, choices=['fasttree', 'iqtree'])
    app.add_argument('--tax_scope_threshold', help='Threshold for assigning the taxonomic scope of a gene cluster [%(default)s]', default=0.9, type=float)
    app.add_argument('--synteny_adj_cutoff', help='Minimum number of times gene clusters should be adjacent to be included in synteny graph [%(default)s]', default=1, type=int)
    app.add_argument('--synteny_mcl_inflation', help='Inflation parameter for mcl clustering of synteny graph. Usually between 1.2 and 5. Higher values produce smaller clusters [%(default)s]', default=5, type=float)

    ########## Workflow ##########
    master = Controller(ap)
    
    wf = Workflow(master)
    
    # If only 1 snakemake command should be run
    if master.snake:
        logging.info('Only running '+master.snake+' snakefile')
        wf.run(snakefile=globals()['WORKFLOW_'+master.snake])
    # If fasta directory is not present
    elif master.fastq:
        # Assembly 
        logging.info('Assemble fastq')
        wf.run(snakefile=WORKFLOW_ASSEMBLY)

    # If reference is present
    elif master.reference: 
        # Quality control
        logging.info('Run quality control for assembly quality and taxonomy verification')
        wf.run(snakefile=WORKFLOW_QUALITY)
    
    elif master.sample_list:
        # Pangenomics
        logging.info('Run pangenome creation and analysis')
        wf.run(snakefile=WORKFLOW_PANGENOME)


if __name__ == '__main__':
    cli()
