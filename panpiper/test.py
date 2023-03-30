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
from Bio import SeqIO
import gzip


class Test(object):
    """
    This class evaluates all user input to prepare for snakemake run.
    """

    def __init__(self, ap):

        args = ap.parse_args()

        for k, v in args.__dict__.items():
            setattr(self, k, v)

        # Logger
        logging.basicConfig(format='\033[36m'+'[%(asctime)s] %(levelname)s:' +
                            '\033[0m'+' %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=self.log_lvl)
        logging.info('Running Panpiper version {}'.format(
            pkg_resources.require("panpiper")[0].version))

        # Force consistency
        self.output = os.path.join(self.output, '')

        # Check cluster info input
        # TODO: update to check cluster information
        if self.cluster:
            if self.profile is None:
                logging.error('A profile folder must be specified')
                sys.exit()
            else: 
                profile = os.path.join(self.profile, "config.yaml")
                if not os.path.isfile(profile):
                    logging.error('config.yaml file not found in profile folder')
                    sys.exit()

        # Check input, params, and output
        self.check_out()
        if not any([self.unlock, self.only_conda, self.snake is not None]):
            if(self.workflow == "assembly"):
                self.check_fastq()
            elif(self.workflow == "quality"):
                self.check_reference()
                self.check_list()
            elif(self.workflow == "pangenome"):
                self.check_list()
                self.check_reference()
            elif(self.workflow == "pyseer"):
                self.check_pyseer()

        self.write_params(args)

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
                if not os.path.join(self.fasta, line, ".fa"):
                    logging.error(
                        'Sample names not present in fasta directory')
                    sys.exit()

    def check_pyseer(self):
        # genes : gene_presence_absence.Rtab check if exists
        if not os.path.isfile(self.genes):
            logging.error('Gene presence/absence file does not exist')
            sys.exit()

        # structure : struct_presence_absence.Rtab check if exists
        if not os.path.isfile(self.structure):
            logging.error('Structure presence/absence file does not exist')
            sys.exit()

        # unitig : unitig.pyseer.gz check that it's zipped
        with gzip.open(self.unitig, 'r') as fh:
            try:
                fh.read(1)
            except gzip.BadGzipFile:
                print('Unitig is a bad gzip file')
                sys.exit()

        # tree : iqtree.nwk check if in newick format
        if not os.path.isfile(self.tree):
            logging.error('Tree file does not exist')
            sys.exit()

        # reference : check if file in text can be found


        # pheno : check if file has sample name as first column

    def write_params(self, args):
        """
        Add parameters to params to snakemake file 
        """

        pars = vars(args)
        self.params = self.output+'parameters.tab'
        fh = open(self.params, 'w')
        for k, v in pars.items():
            fh.write('{}\t{}\n'.format(k, v))
        if not self.reference == "skip":
            ref = os.path.basename(self.reference)
            ref = os.path.splitext(ref)[0]
            fh.write('ref\t{}\n'.format(ref))
        fh.close()
