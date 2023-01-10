import os
import logging
import sys
import pkg_resources
import re
import Path
from Bio import SeqIO

class Controller(object):

    def __init__(self, ap):
     
        args = ap.parse_args()

        for k,v in args.__dict__.items():
            setattr(self, k, v)
        
        # Logger
        logging.basicConfig(format='\033[36m'+'[%(asctime)s] %(levelname)s:'+'\033[0m'+' %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=self.log_lvl)
        logging.info('Running Package version {}'.format(pkg_resources.require("package")[0].version))

        # Force consistency
        self.output = os.path.join(self.output, '')

        # Check cluster info input
        if self.cluster is not None:
            if self.cluster_info is not None:
                if len(re.findall('{cores}|{memory}|{runtime}', self.cluster_info)) != 3 or len(re.findall('{.*?}', self.cluster_info)) != 3: 
                    logging.error('cluster_info has to contain the following special strings: {cores}, {memory}, and {runtime}')
                    sys.exit()
                else:
                    tmp_info = re.sub('{cores}|{memory}|{runtime}','',self.cluster_info)
                    if not bool(re.match('^[a-zA-Z0-9-_ =:,.]+$', tmp_info)):
                        logging.error('Invalid characters in cluster_info')
                        sys.exit()
            else:
                logging.error('cluster_info is required when running on a compute cluster')
                sys.exit()

        # Check input and output
        self.check_params()
        self.check_out()
        if not any([self.unlock, self.only_conda, self.snake is not None]):
            if(self.fastq is not None):
                self.check_fastq()
            if(self.reference is not None):
                if(self.fasta is not None):
                    self.check_reference()
                else:
                    logging.error('Must provide a fasta directory with sample list')
                    sys.exit()
            if(self.sample_list is not None):
                if(self.fasta is not None):
                    self.check_sample()
                else:
                    logging.error('Must provide a fasta directory with sample list')
                    sys.exit()
            else:
                logging.error('Must provide either a fastq directory, reference and fasta directory, or sample list and fasta directory')
                sys.exit()

        self.write_params(args)

    def check_params(self):
        '''
        Return error if parameters are not allowed
        TODO: Add my parameters of interest to this file
        '''
        if self.min_af < 0.5:
            logging.error('min_af lower than 0.5 can lead to unexpected results')
            sys.exit()


    def check_out(self):
        '''
        Creates output directory unless it exists already
        '''
        try:
            os.makedirs(self.output+'logs')
        except FileExistsError:
            logging.warning('Output directory '+self.output+' already exists')

    def check_reference(self):
        '''
        Check the fastq input directory
        '''
        logging.debug('Checking input directory for fastq files')
        cwd = self.fastq.cwd()
        try: 
            for fname in cwd.glob("*.fastq"):
                print('Fastq files found')
        except Exception:
            logging.error('No fastq files found in directory')
            sys.exit()

        # Extract info
        self.fastq_list = set(= [x for x in os.listdir(cwd) if x.endswith(".fastq")])
        

    def check_fasta(self):
        '''
        Check the fasta input directory
        Check the fasta referenc efile
        '''
        logging.debug('Checking input directory for fasta files')
        cwd = self.fastq.cwd()
        try: 
            for fname in cwd.glob("*.fasta"):
                print('Fasta files found')
        except Exception:
            logging.error('No fasta files found in directory')
            sys.exit()

        # Extract info
        self.fasta_list = set(= [x for x in os.listdir(cwd) if x.endswith(".fasta")])
        

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
                elif os.path.join(self.fasta, line, ".fna"):
                elif os.path.join(self.fasta, line, ".fa"):
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
