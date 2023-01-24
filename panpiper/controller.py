# ----------------------------------------------------------------------------
# Copyright (c) 2022--, Panpiper development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import logging
import sys
import pkg_resources
import re
#import Path
from Bio import SeqIO

class Controller(object):
    """
    This class evaluates all user input to prepare for snakemake run.
    """

    def __init__(self, ap):
     
        args = ap.parse_args()

        for k,v in args.__dict__.items():
            setattr(self, k, v)
        
        # Logger
        logging.basicConfig(format='\033[36m'+'[%(asctime)s] %(levelname)s:'+'\033[0m'+' %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=self.log_lvl)
        logging.info('Running Panpiper version {}'.format(pkg_resources.require("panpiper")[0].version))

        # Force consistency
        self.output = os.path.join(self.output, '')

        # Check cluster info input
        # TODO: update to check cluster information
        if self.cluster_type is not None:
            if self.cluster_args is not None:
                if len(re.findall('{cluster.mem}|{cluster.time}|{cluster.cpus}', self.cluster_args)) != 3: 
                    logging.error('cluster_args has to contain the following special strings: {cluster.cpus}, {cluster.mem}, and {cluster.time}')
                    sys.exit()
            else:
                logging.error('cluster_args is required when running on a compute cluster')
                sys.exit()

        # Check input, params, and output
        self.check_params()
        self.check_out()
        if not any([self.unlock, self.only_conda, self.snake is not None]):
            if(self.workflow == "assembly"):
                self.check_fastq()
            elif(self.workflow == "quality"):
                try:
                    self.check_reference()
                    self.check_sample()
                except Exception:
                    logging.error('Must provide reference file and sample list')                    
                    sys.exit()
            elif(self.workflow == "pangenome"):
                try:
                    self.check_sample()
                except Exception:
                    logging.error('Must provide a sample list')
                    sys.exit()

        self.write_params(args)

    def check_params(self):
        """
        Checks parameter validity 

        Raises
        ------
        Error
            'Error: GC content is typically above 20'

        TODO: Add more parameter checks
        """
        if self.gc < 20:
            logging.error('GC content is typically above 20')
            sys.exit()


    def check_out(self):
        """
        Tries to make output directory

        Raises
        ------
        Warning
            'Warning: Output directory already exists'
        """
        try:
            os.makedirs(self.output)
        except FileExistsError:
            logging.warning('Output directory '+self.output+' already exists')

    def check_fastq(self):
        """
        Tries to find input fastq files

        Raises
        ------
        Exception
            'Exception: Fastq files not found'
        """
        logging.debug('Checking input directory for fastq files')
        for file in os.listdir(self.fastq):
            if file.endswith(".fastq"):
                logging.debug('Fastq files found')
                return()
        
        logging.error('No fastq files found in directory')
        sys.exit()        


    def check_reference(self):
        """
        Tries to open reference fasta file 

        Raises
        ------
        Exception
            'Exception: Reference file not found'

        ValueError
            'ValueError: Reference fasta file not in fasta format'
        """
        logging.debug('Checking reference fasta file')

        try: 
            with open(self.reference, "r") as handle:
                fasta = SeqIO.parse(handle, "fasta")
        except Exception:
            logging.error('Reference file not found')
            sys.exit()

        if not (fasta):
            logging.error('Reference fasta file not in fasta format')
            sys.exit()
            
    def check_list(self):
        """
        Checks sample list file for fasta file list 
        """
        logging.debug('Checking sample list file')
        
        with open(self.sample_list, 'r') as fh:
            for line in fh:
                if os.path.join(self.fasta, line, ".fasta"):
                    print('Fasta file found')
                elif os.path.join(self.fasta, line, ".fna"):
                    print('Fasta file found')
                elif os.path.join(self.fasta, line, ".fa"):
                    print('Fasta file found')
                else:
                    logging.error('Sample names not present in fasta directory')
                    sys.exit()

    def write_params(self, args):
        """
        Add parameters to params to snakemake file 
        """
        
        pars = vars(args)
        self.params = self.output+'parameters.tab'
        fh = open(self.params, 'w')
        for k, v in pars.items():
            fh.write('{}\t{}\n'.format(k, v))
        fh.close()
