"""
Authors: Renee Oles
Date Updated: 1/4/2022
Purpose: Quality control
"""
import os
import subprocess
import math
import distutils.util


OUT = config['out']
FASTA = config['fasta']
REFERENCE = config['ref']
PARAMS = config['params']
SAMPLE_LIST = config['list']
SAMPLES_OUT = os.path.join(OUT, 'Quality/sample_list.txt')

with open(PARAMS, 'r') as fh:   
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

def read_complete_files():
    with open(SAMPLE_LIST) as f:
        raw_reads = [sample for sample in f.read().split('\n') if len(sample) > 0]
        return(raw_reads)
READS = read_complete_files()

# Output
rule all:
    input:
        SAMPLES_OUT,


# Remove contigs smaller than 500 bp
# 30 sec per file
rule contig_filter:
    input:
        assembly=os.path.join(FASTA,"{file}/contigs.fa"),
    params:
        contig=500,
    conda:
        "envs/quality.yml"
    log:
        "workflow/report/prokka_filter_{file}.log",
    output:
        fna=os.path.join(OUT,"Assembly_filter/{file}/{file}.fna"),
    shell:
        "reformat.sh in={input.assembly} out={output} minlength={params.contig} &> {log}"


# 5 minutes per sample
rule run_checkm:
    input:
        file=os.path.join(OUT,"Assembly_filter/{file}/{file}.fna"),
    params:
        threads=40,
        binset=os.path.join(OUT,"Assembly_filter/{file}"),
        prefix=os.path.join(OUT,"Quality/Checkm/{file}"),
    output:
        stats=os.path.join(OUT,"Quality/Checkm/{file}/storage/bin_stats.analyze.tsv"),
    conda:
        "envs/quality.yml"
    log:
        log=os.path.join(OUT,"Quality/Checkm/{file}/lineage.log"),
    shell:
        "checkm lineage_wf -t 20 -x fna {params.binset} {params.prefix} &> {log}"


rule checkm_to_graph:
    input:
        stats=expand(os.path.join(OUT,"Quality/Checkm/{file}/storage/bin_stats.analyze.tsv"), file=READS),
    params:
        log=expand(os.path.join(OUT,"Quality/Checkm/{file}/lineage.log"), file=READS),
    output:
        png=os.path.join(OUT,"Quality/checkm_log.txt"),
        stats=os.path.join(OUT,"Quality/checkm_stats.txt"),
    conda:
        "envs/quality.yml"
    log:
        "workflow/report/checkm_graph.log",
    shell:
        """
        Rscript workflow/scripts/checkm-log.R {params.log}  &> {log}
        Rscript workflow/scripts/checkm_bin-stats.R {input.stats}  &> {log}
        """


rule fastani_list_create:
    input:
        fasta=expand(os.path.join(OUT,"Assembly_filter/{file}/{file}.fna"), file=READS),
    params:
        ref=REFERENCE,
    output:
        os.path.join(OUT,"Quality/fastani_list.txt"),
    log:
        "workflow/report/fastani_list_create.log",
    shell:
        """
        chmod u+x workflow/scripts/create_list.sh
        workflow/scripts/create_list.sh {input} {params} &> {log}
        """


rule fastani:
    input:
        os.path.join(OUT,"Quality/fastani_list.txt"),
    output:
        os.path.join(OUT,"Quality/matrix.txt"),
    conda:
        "envs/quality.yml"
    log:
        "workflow/report/fastani.log",
    shell:
        "fastANI --ql {input} --rl {input} -o {output} -t 20 &> {log}"


rule fastani_long_to_wide:
    input:
        os.path.join(OUT,"Quality/matrix.txt"),
    output:
        txt=os.path.join(OUT,"Quality/matrix_wide.txt"),
    log:
        "workflow/report/fastani_long_to_wide.log",
    shell:
        "python workflow/scripts/long_to_wide.py {input} {output.txt} &> {log}"


# Filter files based off user-defined rules
rule filter_files:
    input:
        ani=os.path.join(OUT,"Quality/matrix_wide.txt"),
        stat=os.path.join(OUT,"Quality/checkm_stats.txt"),
    params:
        checkm=os.path.join(OUT,"Quality/checkm_log.txt"),
        ref=REFERENCE,
        gc=param_dict["gc"],
        genome_size=param_dict["genome_size"],
    log:
        "workflow/report/filter_files.log",
    output:
        os.path.join(OUT,"Quality/sample_list.txt"),
    shell:
        """
        python workflow/scripts/filter_isolates.py -o {output} -a {input.ani} -l {params.checkm} -s {input.stat} -r {params.ref} -gc {params.gc} -g {params.genome_size} &> {log}
        """