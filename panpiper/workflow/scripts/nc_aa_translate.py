# ----------------------------------------------------------------------------
# Copyright (c) 2022--, Panpiper development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os

def translate(fasta, output):
    with open(output, 'w') as aa_fa:
        for dna_record in SeqIO.parse(fasta, 'fasta'):
            dna_seqs = [dna_record.seq]
            # generate all translation frames
            aa_seqs = (s[i:].translate(to_stop=True) for i in range(3) for s in dna_seqs)
            # select the longest one
            max_aa = max(aa_seqs, key=len)
            # write new record
            aa_record = SeqRecord(max_aa, id=dna_record.id, description="")
            SeqIO.write(aa_record, aa_fa, 'fasta')

if __name__ == "__main__":
    fasta = sys.argv[1]
    output = sys.argv[2]
    translate(fasta, output)