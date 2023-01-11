import sys
import yaml
import re
import os
import subprocess
import logging
import glob

class Message(object):

    def __init__(self, obj):
        for key, val in vars(obj).items():
            setattr(self, key, val)
        
        self.commands = []

    def add(self, x):
        
        self.commands.append(x)

    def check(self):
        
        # Get latest command
        last = self.commands[-1]

        if 'workflow/assembly.smk' in last:
            logging.info('Finished assembly')

        if 'workflow/quality.smk' in last:
            logging.debug('Filtering done')

        if 'workflow/pangenome.smk' in last:
            logging.debug('Pangenome created')