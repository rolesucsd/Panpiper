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