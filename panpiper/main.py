#!/usr/bin/env python
# ----------------------------------------------------------------------------
# Copyright (c) 2022--, Panpiper development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

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

from panpiper.test import Test
from panpiper.workflow import Workflow

# Snakefiles
_ROOT = os.path.abspath(os.path.dirname(__file__))
WORKFLOW_ASSEMBLY = os.path.join(_ROOT, 'workflow', 'assembly.smk')
WORKFLOW_QUALITY = os.path.join(_ROOT, 'workflow', 'quality.smk')
WORKFLOW_PANGENOME = os.path.join(_ROOT, 'workflow', 'pangenome.smk')
WORKFLOW_PYSEER = os.path.join(_ROOT, 'workflow', 'pyseer.smk')


def cli():

    ########## Arguments ##########
    ap = argparse.ArgumentParser(description='Package version ', add_help=False,
                                 allow_abbrev=False, formatter_class=RawTextHelpFormatter)

    # Required
    apr = ap.add_argument_group('main arguments')
    apr.add_argument('-q', '--fastq', help='Path to fastq file directory',
                     required=False, default="skip", type=str)
    apr.add_argument('-a', '--fasta', help='Path to fasta file directory',
                     required=False, default="skip", type=str)
    apr.add_argument('-s', '--sample_list', help='Line delimited file of samples passing quality control',
                     required=False, default="skip", type=str)
    apr.add_argument('-g', '--genes', help='File of gene presence/absence output from Panaroo/Roary',
                     required=False, default="skip", type=str)
    apr.add_argument('-p', '--structure', help='File of structure presence/absence output from Panaroo/Roary',
                     required=False, default="skip", type=str)
    apr.add_argument('-u', '--unitig', help='File of unitig list',
                     required=False, default="skip", type=str)
    apr.add_argument('-t', '--tree', help='Newick tree of all samples from IQTree',
                     required=False, default="skip", type=str)
    apr.add_argument('-o', '--output',
                     help='Prefix for output directory', required=True)
    apr.add_argument('-r', '--reference', help='Reference fasta file, includes path',
                     required=False, default="skip", type=str)
    apr.add_argument('-w', '--workflow', help='Choose which workflow to run [%(default)s]',
                     default='assembly', type=str, choices=['assembly', 'quality', 'pangenome', 'pyseer'])
    apr.add_argument(
        '-f', '--pheno', help='The filename and full path of the phenotype file', default="skip", type=str)

    # Cluster arguments
    apc = ap.add_argument_group('compute cluster arguments')
    apc.add_argument(
        '--cluster_type', help='Cluster compute structure [%(default)s]', default=None, type=str, choices=[None, 'qsub', 'slurm'])
    apc.add_argument('--cluster_config',
                     help='Cluster json file [%(default)s]', default=None, type=str)
    apc.add_argument(
        '--cluster_args', help='Cluster scheduler arguments when submitting cluster jobs.\nFollow the format: "sbatch -A {cluster.account} --mem {cluster.mem} -t {cluster.time} --cpus-per-task {cluster.cpus}"', default=None, type=str)
    apc.add_argument(
        '--max_jobs', help='Maximum number of cluster jobs [%(default)s]', default=40, type=int)

    # Optional
    apo = ap.add_argument_group('optional arguments')
    apo.add_argument("-h", "--help", action="help",
                     help="show this help message and exit")
    apo.add_argument(
        "-j", "--jobs", help="Number of jobs to run at a time", default="50", type=str)
    apo.add_argument(
        '--max_cores', help='Maximum number of cores [%(default)s]', default=20, type=int)
    apo.add_argument(
        '--max_mem', help='Maximum memory in GB [%(default)s]', default=250, type=int)
    apo.add_argument('--log_lvl', help='Logging level [%(default)s].',
                     default='INFO', type=str, choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'])
    apo.add_argument(
        '--only_conda', help='Only install conda environments, then exit', action='store_true')
    apo.add_argument(
        '--snake', help='Only run specific snakemake command. For debug purposes', type=str)
    apo.add_argument(
        '--dry_run', help='Only run specific snakemake command. For debug purposes', action='store_true')
    apo.add_argument(
        '--unlock', help='Unlock snakemake directory in case of unexpected exists, then exit', action='store_true')

    # Parameters
    # Defaults set for Bacteroides fragilis
    app = ap.add_argument_group('parameters')
    app.add_argument(
        '--gc', help='GC content of species of interest', default=43.19, type=float)
    app.add_argument(
        '--genome_size', help='Genome size of species of interest', default=5205140, type=int)
    app.add_argument('--ref', help='Name of reference',
                     default="9343", type=str)
    app.add_argument(
        '--ani_cutoff', help='Average percent identity to reference cutoff', default=95, type=float)
    app.add_argument('--contig_number',
                     help='Max number of contigs', default=1000, type=int)
    app.add_argument('--n50', help='N50 cutoff', default=5000, type=int)
    app.add_argument('--eggnog_dir', help='Path to the eggnog database directory',
                     default="panpiper/databases/eggnog", type=str)
    app.add_argument('--kraken_dir', help='Path to the kraken2 database directory',
                     default="panpiper/databases/kraken", type=str)
    app.add_argument('--bakta_dir', help='Path to the bakta database directory',
                     default="panpiper/databases/bakta", type=str)
    app.add_argument(
        '--pheno_column', help='The column in the phenotype file to use for the association study', type=str)

    ########## Workflow ##########
    master = Test(ap)

    wf = Workflow(master)

    # If only 1 snakemake command should be run
    if master.snake:
        logging.info('Only running '+master.snake+' snakefile')
        wf.run(snakefile=globals()['WORKFLOW_'+master.snake])

    # If workflow is set to Assembly
    elif (master.workflow == "assembly"):
        # Assembly
        logging.info('Assemble fastq')
        wf.run(snakefile=WORKFLOW_ASSEMBLY)

    # If workflow is set to Quality
    elif (master.workflow == "quality"):
        # Quality control
        logging.info(
            'Run quality control for assembly quality and taxonomy verification')
        logging.info('The parameters used in this study are- GC:' + str(master.gc) + ', Genome Size:' + str(master.genome_size) + ', Reference:' + master.ref +
                     ', ANI Cutoff:' + str(master.ani_cutoff) + ', Contig Number:' + str(master.contig_number) + ', N50:' + str(master.n50))
        wf.run(snakefile=WORKFLOW_QUALITY)

    # If workflow is set to Pangenome
    elif (master.workflow == "pangenome"):
        logging.info('Run pangenome creation and analysis')
        logging.info('The parameters used in this study are- Eggnog Database Directory: ' + master.eggnog_dir +
                     ', Kraken Database Directory: ' + master.kraken_dir + ', Bakta Database Directory: ' + master.bakta_dir)
        wf.run(snakefile=WORKFLOW_PANGENOME)

    # If workflow is set to Pyseer
    elif (master.workflow == "pyseer"):
        logging.info('Run genome-wide association study')
        logging.info(
            'The parameters used in this study are- Phenotype column:' + str(master.pheno_column))
        wf.run(snakefile=WORKFLOW_PYSEER)


if __name__ == '__main__':
    cli()
