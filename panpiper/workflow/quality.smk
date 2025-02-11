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
SAMPLES_OUT = os.path.join(OUT, 'Quality/quality_report.html')
PATH = config['scripts']

with open(PARAMS, 'r') as fh:   
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

PARAMS_REF = param_dict["ref"]
ANI_CUTOFF = param_dict["ani_cutoff"]
CONTIG_NUMBER = param_dict["contig_number"]
N50 = param_dict["n50"]

def read_complete_files():
    with open(SAMPLE_LIST) as f:
        raw_reads = [sample for sample in f.read().split('\n') if len(sample) > 0]
        return(raw_reads)
READS = read_complete_files()

def find_assembly(wildcards):
    for ext in ['fa', 'fna', 'fasta']:
        path = os.path.join(FASTA, f"{wildcards.file}.{ext}")
        if os.path.isfile(path):
            return path
    raise ValueError(f"No assembly file found for {path} with extensions .fa, .fna, or .fasta")

# Output
rule all:
    input:
        SAMPLES_OUT,

rule contig_filter:
    input:
        assembly=find_assembly,
    params:
        contig=500,
    conda:
        "envs/bbmap.yml"
    log:
        os.path.join(OUT,"report/prokka_filter_{file}.log"),
    resources:
        mem="1G"
    threads: 1
    benchmark:
        os.path.join(OUT,"benchmark/prokka_filter_{file}.benchmark"),    
    output:
        fna=os.path.join(OUT,"Quality/Samples/{file}/{file}.fna"),
    shell:
        "reformat.sh in={input.assembly} out={output} minlength={params.contig} &> {log}"

rule run_checkm2:
    input:
        file=os.path.join(OUT,"Quality/Samples/{file}/{file}.fna"),
    params:
        binset=os.path.join(OUT,"Quality/Samples/{file}/checkm"),
        checkmdb="/ddn_scratch/roles/Panpiper/panpiper/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd"
    conda:
        "envs/checkm2.yml"
    log:
        log=os.path.join(OUT,"report/checkm2_{file}.log"),
    benchmark:
        os.path.join(OUT,"benchmark/checkm2_{file}.benchmark"),
    output:
        stats=os.path.join(OUT,"Quality/Samples/{file}/checkm/quality_report.tsv"),
    shell:
        """
        checkm2 predict --threads 20 --input {input} --database_path {params.checkmdb} --output-directory {params.binset} --force &> {log}
        """

rule checkm_to_graph:
    input:
        stats=expand(os.path.join(OUT,"Quality/Samples/{file}/checkm/quality_report.tsv"), file=READS),
    params:
        outpref=OUT,
        path=PATH,
    log:
        os.path.join(OUT,"report/checkm_cat.log"),
    benchmark:
        os.path.join(OUT,"benchmark/checkm_graph.benchmark"),
    resources:
        mem="1G"
    threads: 1        
    output:
        stats=os.path.join(OUT,"Quality/CheckM/checkm_cat.txt"),
    shell:
        """
        python {params.path}/checkm_cat.py {params.outpref} {input.stats} &> {log}
        """

rule fastani:
    input:
        ref=REFERENCE,
        file=os.path.join(OUT,"Quality/Samples/{file}/{file}.fna"),
    conda:
        "envs/fastani.yml"
    log:
        os.path.join(OUT,"report/fastani_{file}.log"),
    resources:
        mem="1G"
    threads: 1
    benchmark:
        os.path.join(OUT,"benchmark/fastani_{file}.benchmark"),    
    output:
        os.path.join(OUT,"Quality/FastANI/{file}.txt"),
    shell:
        "fastANI -q {input.file} -r {input.ref} -o {output} &> {log}"


rule fastani_concat:
    input:
        expand(os.path.join(OUT,"Quality/FastANI/{file}.txt"), file=READS),
    resources:
        mem="1G"
    threads: 1
    output:
        txt=os.path.join(OUT,"Quality/FastANI/fastani_summary.txt"),
    shell:
        "cat {input} > {output}"


# Filter files based off user-defined rules
# Need to edit this python file bc it's messy rn
rule filter_files:
    input:
        ani=os.path.join(OUT,"Quality/FastANI/fastani_summary.txt"),
        checkm=os.path.join(OUT,"Quality/CheckM/checkm_cat.txt"),
    params:
        outpath=os.path.join(OUT,"Quality"),
        ref=PARAMS_REF,
        ac=ANI_CUTOFF,
        n=N50,
        path=PATH,
    log:
        os.path.join(OUT,"report/filter_files.log"),
    benchmark:
        os.path.join(OUT,"benchmark/filter_files.benchmark"),    
    resources:
        mem="1G"
    threads: 1
    output:
        os.path.join(OUT,"Quality/sample_list.txt"),
    shell:
        """
        python {params.path}/filter_isolates.py -o {params.outpath} -a {input.ani} -k {input.checkm} -r {params.ref} -ac {params.ac} -n {params.n} &> {log}
        """

rule print_results:
    input:
        ani=os.path.join(OUT,"Quality/FastANI/fastani_summary.txt"),
        checkm=os.path.join(OUT,"Quality/CheckM/checkm_cat.txt"),
        passed=os.path.join(OUT,"Quality/sample_list.txt"),
    params:
        ref=PARAMS_REF,
        outdir=os.path.join(OUT, 'Quality')
    conda:
        "envs/r.yml"
    log:
        os.path.join(OUT,"report/print_results.log"),
    benchmark:
        os.path.join(OUT,"benchmark/print_results.benchmark"),    
    resources:
        mem="5G"
    threads: 1
    output:
        SAMPLES_OUT,
    shell:
        "Rscript -e \"rmarkdown::render('panpiper/workflow/scripts/quality_report.Rmd', output_dir = '{params.outdir}', params=list(checkm = '{input.checkm}', ani = '{input.ani}', passed = '{input.passed}', ref = '{params.ref}'))\" &> {log}"


