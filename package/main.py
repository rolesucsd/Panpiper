#!/usr/bin/env python

import argparse
import logging
import os
import sys
import pkg_resources
# Delete after package creation
import importlib

from os.path import join, dirname, basename
from shutil import copyfile
from argparse import RawTextHelpFormatter

# Add back after package creation
#from package.controller import Controller
#from package.workflow import Workflow

# Delete after package creation
from controller import Controller
from workflow import Workflow


# Snakefiles
_ROOT = os.path.abspath(os.path.dirname(__file__))
WORKFLOW_ASSEMBLY = os.path.join(_ROOT, 'workflow', 'assembly.smk')
WORKFLOW_QUALITY = os.path.join(_ROOT, 'workflow', 'quality.smk')
WORKFLOW_PANGENOME = os.path.join(_ROOT, 'workflow', 'pangenome.smk')

def cli():
    
    ########## Arguments ##########
#    ap = argparse.ArgumentParser(description='Package version {}'.format(pkg_resources.require("package")[0].version), add_help=False, allow_abbrev=False, formatter_class=RawTextHelpFormatter)
    ap = argparse.ArgumentParser(description='Package version ', add_help=False, allow_abbrev=False, formatter_class=RawTextHelpFormatter)
    
    # Required
    apr = ap.add_argument_group('main arguments')
    apr.add_argument('-q', '--fastq', help='Path to fastq file directory', required=False, default="skip", type=str)
    apr.add_argument('-a', '--fasta', help='Path to fasta file directory', required=False, default="skip", type=str)
    apr.add_argument('-s', '--sample_list', help='Line delimited file of samples passing quality control', required=False, default="skip", type=str)
    apr.add_argument('-o', '--output', help='Prefix for output directory', required=True)
    apr.add_argument('-r', '--reference', help='Reference fasta file, includes path', required=False, default="skip", type=str)

    # Cluster arguments
    apc = ap.add_argument_group('compute cluster arguments')
    apc.add_argument('--cluster', help='Cluster compute structure [%(default)s]', default=None, type=str, choices=[None,'qsub','slurm','drmaa'])
    apc.add_argument('--cluster_info', help='Cluster scheduler arguments when submitting cluster jobs.\nHas to contain the following special strings:\n{memory}, {cores}, and {runtime}.\nThese special strings will be substituted by maginator to indicate resources for each job.\n{memory} is substituted for the memory in GB.\n{runtime} is substituted with the time in the following format: DD:HH:MM:SS.\nCan also contain user names, groups, etc. required by the cluster scheduler', default=None, type=str)
    apc.add_argument('--max_jobs', help='Maximum number of cluster jobs [%(default)s]', default=50, type=int)
    
    # Optional
    apo = ap.add_argument_group('optional arguments')
    apo.add_argument("-h", "--help", action="help", help="show this help message and exit")
    apo.add_argument("-j", "--jobs", help="Number of jobs to run at a time", default="1", type=str)
    apo.add_argument('--max_cores', help='Maximum number of cores [%(default)s]', default=40, type=int)
    apo.add_argument('--max_mem', help='Maximum memory in GB [%(default)s]', default=180, type=int)
    apo.add_argument('--log_lvl', help='Logging level [%(default)s].', default='INFO', type=str, choices=['DEBUG','INFO','WARNING','ERROR'])
    apo.add_argument('--only_conda', help='Only install conda environments, then exit', action='store_true')
    apo.add_argument('--snake', help='Only run specific snakemake command. For debug purposes', type=str)
    apo.add_argument('--unlock', help='Unlock snakemake directory in case of unexpected exists, then exit', action='store_true')

    # Parameters
    ## Defaults set for Bacteroides fragilis
    app = ap.add_argument_group('parameters')
    app.add_argument('--gc', help='GC content of species of interest', default=43.19, type=float)
    app.add_argument('--genome_size', help='Genome size of species of interest', default = 5205140, type=int)
    app.add_argument('--ref', help='Name of reference', default = "9343", type=str)

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
