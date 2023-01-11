import sys
import yaml
import re
import os
import subprocess
import logging
import glob

from messages import Message

class Workflow(object):

    def __init__(self, obj):
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run(self, snakefile):

        # Start message instance
        messages = Message(self)

        # Define core snakemake command
        # need to edit j
        cmd = ['snakemake',
               '--use-conda',
               '--latency-wait', '20',
               '-s', snakefile,
               '--config',
               'out='+self.output,
               'fastq='+self.fastq,
               'fasta='+self.fasta,
               'list='+self.sample_list,
               'ref='+self.reference,
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
                    os.mkdir(self.output + 'report/cluster_err')
                except FileExistsError:
                    pass
                try:
                    os.mkdir(self.output + 'report/cluster_out')
                except FileExistsError:
                    pass
            else:
                try:
                    os.mkdir(self.output + 'report/drmaa')
                except FileExistsError:
                    pass
            
            # Add cluster info to snakemake command
            cmd += ['--jobs', str(self.max_jobs),
                    '--local-cores', str(self.max_cores)]
                
        if self.cluster_type == 'qsub':
            
            cluster_cmd = '--cluster-config' + ' ' + self.cluster_config
#            cluster_cmd = cluster_cmd + ' -e ' + self.output + 'report/cluster_err' + ' -o ' + self.output + 'report/cluster_out'
            cluster_cmd = cluster_cmd + ' --cluster ' + '"' + self.cluster_args + '"'
            # Final snakemake command
            cmd += [cluster_cmd]

        if self.cluster_type == 'slurm':
            
            cluster_cmd = '--cluster-config' + ' ' + self.cluster_config
#            cluster_cmd = cluster_cmd + ' -e ' + self.output + 'report/cluster_err' + ' -o ' + self.output + 'report/cluster_out'
            cluster_cmd = cluster_cmd + ' --cluster ' + '"' + self.cluster_args + '"'
            # Final snakemake command
            cmd += [cluster_cmd]

        # Only install conda envs if only_conda
        if self.only_conda:
            logging.info('Only creating conda environments.')
            cmd.append(' --conda-create-envs-only')

        # If unlocking
        if self.unlock:
            logging.info('Unlocking working directory.')
            cmd.append('--unlock')

        logging.info(' '.join(cmd))
        # Start snakemake process and read stdout and stderr (also save in logger)
        process = subprocess.Popen(cmd, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, text=True)
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

        # Save information on what has been run
        messages.add(snakefile) 

        # If unlocking exit program
        if self.unlock:
            sys.exit()

        # Print info
        if not self.only_conda:
            messages.check()
