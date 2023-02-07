#!/usr/bin/env python
# ----------------------------------------------------------------------------
# Copyright (c) 2023--, Panpiper development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
from io import StringIO
import json
from unittest import TestCase

import numpy as np
import numpy.testing as npt
import pytest

from panpiper.main import Main

__author__ = "Renee Oles"
__copyright__ = "Copyright 2023, Panpiper development team"
__credits__ = ["Renee Oles"]
__license__ = "BSD"
__maintainer__ = "Renee Oles"
__email__ = "roles@health.ucsd.edu"

class MainTest(TestCase):

    """Tests of parse functions"""
    def setUp(self):
        """define some top-level data"""
        self.outpath = outpath # set to test_data dir 
        self.fastq = fastq # test_data/fastq
        self.assembly1 = assembly1 # test_data/fasta1
        self.assembly2 = assembly2 # test_data/fasta2
        self.list1 = list1 # test_data/files.txt
        self.list2 = list2 # test_data/sample_list.txt
        self.reference = reference # test_data/9343.fna
        self.reference_empty = reference_empty # empty fasta file

        def test_unlock(self):
            exp = 0
            obs = "panpiper -o " +self.outpath+ " --unlock" 
            self.assertEqual(obs, exp)

        def test_assembly(self):
            exp = 0
            obs = "panpiper --dry_run -o " +self.outpath+ " -w assembly -q " +self.fastq
            self.assertEqual(obs, exp)

        def test_quality(self):
            exp = 0
            obs = "panpiper --dry_run -o " + self.outpath+ " -a " +self.assembly1+ " -s "  +self.list1+ " -r " +self.reference+ "-w quality" 
            self.assertEqual(obs, exp)

        def test_pangenome(self):
            exp = 0
            obs = "panpiper --dry_run -o " + self.outpath+ " -a " +self.assembly2+ " -s "  +self.list2+ " -r " +self.reference+ "-w pangenome"
            self.assertEqual(obs, exp)

        def test_wrong_assembly_quality(self):
            exp = 0
            obs = "panpiper --dry_run -o " + self.outpath+ " -a " +self.assembly2+ " -s " +self.list1+ " -r " +self.reference+ "-w quality" 
            self.assertEqual(obs, exp)

        def test_wrong_assembly_pangenome(self):
            exp = 0
            obs = "panpiper --dry_run -o " + self.outpath+ " -a " +self.assembly1+ " -s " +self.list2+ " -r " +self.reference+ "-w pangenome"
            self.assertEqual(obs, exp)

        def test_no_outpath(self):
            exp = 0
            obs = "panpiper -w assembly -q" +self.fastq  
            self.assertEqual(obs, exp)

        def test_incorrect_input(self):
            exp = 0
            obs = "panpiper --dry_run -o " +self.outpath+ " -q " +self.fastq+ " -s "  +self.list1+ " -r " +self.reference+ " -w quality" 
            self.assertEqual(obs, exp)

        def test_no_refernce(self):
            exp = 0
            obs = "panpiper --dry_run -o " +self.outpath+ " -a " +self.assembly1+ " -s " +self.list1+ " -w quality"
            self.assertEqual(obs, exp)

        def test_empty_refernce(self):
            exp = 0
            obs = "panpiper --dry_run -o " +self.outpath+ " -q " +self.fastq+ " -s "  +self.list1+ " -r " +self.reference_empty+ " -w quality"
            self.assertEqual(obs, exp)

outpath = "test_data"
fastq = "test_data/fastq"
assembly1 = "test_data/fasta1"
assembly2 = "test_data/fastaq"
list1 = "test_data/files.txt"
list2 = "test_data/sample_list.txt"
reference = "test_data/9343.fna"
reference_empty = "test_data/empty.fna"

if __name__ == '__main__':
    main()

