import os
import logging
import sys
import pkg_resources
import re
#import Path
from Bio import SeqIO

class Controller(object):

    def __init__(self, ap):
     
        args = ap.parse_args()

        for k,v in args.__dict__.items():
            setattr(self, k, v)
        
        # Logger
        logging.basicConfig(format='\033[36m'+'[%(asctime)s] %(levelname)s:'+'\033[0m'+' %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=self.log_lvl)
#        logging.info('Running Package version {}'.format(pkg_resources.require("package")[0].version))

        # Force consistency
        self.output = os.path.join(self.output, '')

        # Check cluster info input
        if self.cluster_type is not None:
            if self.cluster_args is not None:
                if len(re.findall('{cluster.mem}|{cluster.time}|{cluster.cpus}', self.cluster_args)) != 3: 
                    logging.error('cluster_args has to contain the following special strings: {cluster.cpus}, {cluster.mem}, and {cluster.time}')
                    sys.exit()
            else:
                logging.error('cluster_args is required when running on a compute cluster')
                sys.exit()

        # Check input and output
        self.check_params()
        self.check_out()
        if not any([self.unlock, self.only_conda, self.snake is not None]):
            if(self.workflow is "assembly"):
                self.check_fastq()
            elif(self.workflow is "quality"):
                try:
                    self.check_reference()
                    self.check_sample()
                except Exception:
                    logging.error('Must provide reference file and sample list')                    
                    sys.exit()
            elif(self.workflow is "pangenome"):
                try:
                    self.check_sample()
                except Exception:
                    logging.error('Must provide a sample list')
                    sys.exit()

        self.write_params(args)

    def check_params(self):
        '''
        Return error if parameters are not allowed
        TODO: Add my parameters of interest to this file
        '''
        if self.gc < 20:
            logging.error('GC content is typically above 20')
            sys.exit()


    def check_out(self):
        '''
        Creates output directory unless it exists already
        '''
        try:
            os.makedirs(self.output+'logs')
        except FileExistsError:
            logging.warning('Output directory '+self.output+' already exists')

    def check_fastq(self):
        '''
        Check the fastq input directory
        '''
        logging.debug('Checking input directory for fastq files')
        for file in os.listdir(self.fastq):
            if file.endswith(".fastq"):
                logging.debug('Fastq files found')
                return()
        
        logging.error('No fastq files found in directory')
        sys.exit()        


    def check_reference(self):
        '''
        Check the fasta reference file
        '''
        logging.debug('Checking reference fasta file')

        try: 
            with open(self.reference, "r") as handle:
                fasta = SeqIO.parse(handle, "fasta")
        except Exception:
            logging.error('Reference file not found')
            sys.exit()

        if not (fasta):
            logging.error('Reference fasta file is not correct format')
            sys.exit()
            
    def check_list(self):
        '''
        Read sample list file
        '''
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
        '''
        Write parameters for snakemake workflows to a file
        '''

        pars = vars(args)
        self.params = self.output+'parameters.tab'
        fh = open(self.params, 'w')
        for k, v in pars.items():
            fh.write('{}\t{}\n'.format(k, v))
        fh.close()
