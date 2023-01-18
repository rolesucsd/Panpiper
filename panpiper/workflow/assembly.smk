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
PARAMS = config['params']
ASSEMBLY_OUT = os.path.join(OUT, 'complete_assembly_files.txt')
(READS,) = glob_wildcards(os.path.join(FASTQ,"{file}_2.fastq"))

with open(PARAMS, 'r') as fh:   
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

rule all:
    input:
        ASSEMBLY_OUT,


# Assembly with shovill https://github.com/tseemann/shovill
## Includes trimming and error correction
## nocorr just doesn't do pilon, it will still filter low coverage or small contig size --nocorr
rule shovill:
    input:
        read1=os.path.join(FASTQ, '{file}_1.fastq'),
        read2=os.path.join(FASTQ, '{file}_2.fastq'),
    params:
        prefix=os.path.join(OUT, 'Assembly/{file}'),
        tmp="tmp",
    log:
        os.path.join(OUT,"report/shovill_{file}.log"),
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
        os.path.join(OUT, 'Assembly/')
    log:
        os.path.join(OUT,"report/assembly_complete.log"),
    output:
        ASSEMBLY_OUT,
    shell:
        """
        chmod u+x workflow/scripts/move_files.sh
        workflow/scripts/move_files.sh {params} &> {log}
        """
