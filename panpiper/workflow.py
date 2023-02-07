# ----------------------------------------------------------------------------
# Copyright (c) 2022--, Panpiper development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import yaml
import re
import os
import subprocess
import logging
import glob
import shlex
import subprocess


class Workflow(object):

    def __init__(self, obj):
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run(self, snakefile):
        """
        Construct the snakemake command based off user-defined parameters
        Will produce a command, cmd, which will be launched

        Parameters 
        ----------
        self: object, required
        Contains the parameters defined by the user

        snakefile: string, required
        The snakefile selected for the desired operation
        """
        # Define core snakemake command
        # need to edit j
        cmd = ['snakemake',
               '--use-conda',
               '--latency-wait', '20',
               '--rerun-incomplete',
               '-s', snakefile,
               '--config',
               'out='+self.output,
               'fastq='+self.fastq,
               'fasta='+self.fasta,
               'list='+self.sample_list,
               'ref='+self.reference,
               'structure='+self.structure,
               'unitig='+self.unitig,
               'tree='+self.tree,
               'genes='+self.genes,
               'pheno='+self.pheno,
               'params='+self.params]

        # If run on server
        if self.cluster_config == None:
            cmd += ['--cores', str(self.max_cores)]

        # If run on a cluster
        else:

            # Make logging dirs
            if self.cluster_type in ('qsub', 'slurm'):
                try:
                    os.mkdir(self.output + 'report')
                except FileExistsError:
                    pass

            # Add cluster info to snakemake command
            cmd += ['--jobs', str(self.max_jobs),
                    '--local-cores', str(self.max_cores)]

        if self.cluster_type == 'qsub' or self.cluster_type == 'slurm':

            cluster_cmd = ['--cluster-config', self.cluster_config]
            cluster_args_mod = '"' + self.cluster_args + '"'
            cluster_cmd += [' --cluster ', cluster_args_mod]
            # Final snakemake command
            cmd += cluster_cmd

        # Only install conda envs if only_conda
        if self.only_conda:
            logging.info('Only creating conda environments.')
            cmd.append(' --conda-create-envs-only')

        if self.dry_run:
            logging.info('Creating a dry run.')
            cmd.append('-n')
        # If unlocking

        if self.unlock:
            logging.info('Unlocking working directory.')
            cmd.append('--unlock')

        logging.info(' '.join(cmd))
        args = shlex.split(' '.join(cmd))

        # Run the snakemake workflow
        process = subprocess.Popen(
            args, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, text=True)
        if self.log_lvl == 'DEBUG':
            for line in iter(process.stdout.readline, ""):
                logging.debug(line.strip())

        process.wait()
        logging.debug('Snakemake returncode: '+str(process.returncode))

        # Check returncode
        if process.returncode != 0:
            for line in iter(process.stdout.readline, ""):
                logging.error(line.strip())
            logging.error('Snakemake returncode: '+str(process.returncode))
            sys.exit()

        # If unlocking exit program
        if self.unlock:
            sys.exit()
