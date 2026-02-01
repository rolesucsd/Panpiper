# ----------------------------------------------------------------------------
# Copyright (c) 2022--, Panpiper development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""
Authors: Renee Oles
Date Updated: 1/4/2022
Purpose: Assembly
"""
import os
import subprocess
import math
import distutils.util


OUT = config['out']
FASTQ = config['fastq']
PATH = config['scripts']
ASSEMBLY_OUT = os.path.join(OUT, 'Assembly/complete_assembly_files.txt')
(READS,) = glob_wildcards(os.path.join(FASTQ,"{file}_2.fastq.gz"))

rule all:
    input:
        ASSEMBLY_OUT,


# Assembly with shovill https://github.com/tseemann/shovill
# Time: 5-40 minutes per sample
rule shovill:
    input:
        read1=os.path.join(FASTQ, '{file}_1.fastq.gz'),
        read2=os.path.join(FASTQ, '{file}_2.fastq.gz'),
    params:
        prefix=os.path.join(OUT, 'Assembly/{file}'),
        tmp="tmp",
    log:
        os.path.join(OUT,"report/shovill_{file}.log"),
    benchmark:
        os.path.join(OUT,"benchmark/shovill_{file}.benchmark"),
    conda:
        "envs/shovill.yml"
    output:
        os.path.join(OUT, 'Assembly/{file}/contigs.fa'),
    shell:
        """
        shovill --nocorr --trim --force --cpus 0 --ram 8 --opts '-m 1024' --outdir {params.prefix} --R1 {input.read1} --R2 {input.read2} &> {log} || touch {output}
        """

# Moves files that were successfully assembled to the input folder for use in the next workflow
rule move_complete:
    input:
        expand(os.path.join(OUT, 'Assembly/{file}/contigs.fa'), file=READS),
    params:
        outdir=os.path.join(OUT, 'Assembly/'),
        script=os.path.join(PATH, 'move_files.sh'),
    log:
        os.path.join(OUT,"report/move_complete.log"),
    benchmark:
        os.path.join(OUT,"benchmark/move_complete.benchmark"),
    output:
        ASSEMBLY_OUT,
    shell:
        """
        chmod u+x {params.script}
        {params.script} {params.outdir} &> {log}
        """
