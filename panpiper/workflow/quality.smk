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

# Output
rule all:
    input:
        SAMPLES_OUT,


# Remove contigs smaller than 500 bp
# 30 sec per file
rule contig_filter:
    input:
        assembly=os.path.join(FASTA,"{file}.fa"),
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
        fna=os.path.join(OUT,"Quality/Assembly_filter/{file}/{file}.fna"),
    shell:
        "reformat.sh in={input.assembly} out={output} minlength={params.contig} &> {log}"


# 5 minutes per sample
rule run_checkm:
    input:
        file=os.path.join(OUT,"Quality/Assembly_filter/{file}/{file}.fna"),
    params:
        binset=os.path.join(OUT,"Quality/Assembly_filter/{file}"),
    conda:
        "envs/checkm.yml"
    log:
        log=os.path.join(OUT,"Quality/Assembly_filter/{file}/lineage.log"),
    benchmark:
        os.path.join(OUT,"benchmark/checkm_{file}.benchmark"),    
    resources:
        mem="50G"
    threads: 20
    output:
        stats=os.path.join(OUT,"Quality/Assembly_filter/{file}/storage/bin_stats.analyze.tsv"),
    shell:
        """
        checkm lineage_wf -t 20 -x fna {params.binset} {params.binset} &> {log}
        """


rule checkm_to_graph:
    input:
        stats=expand(os.path.join(OUT,"Quality/Assembly_filter/{file}/storage/bin_stats.analyze.tsv"), file=READS),
    params:
        log=expand(os.path.join(OUT,"Quality/Assembly_filter/{file}/lineage.log"), file=READS),
    conda:
        "envs/r.yml"
    log:
        os.path.join(OUT,"report/checkm_graph.log"),
    benchmark:
        os.path.join(OUT,"benchmark/checkm_graph.benchmark"),
    resources:
        mem="1G"
    threads: 1        
    output:
        png=os.path.join(OUT,"Quality/CheckM/checkm_log.txt"),
        stats=os.path.join(OUT,"Quality/CheckM/checkm_stats.txt"),
    shell:
        """
        Rscript panpiper/workflow/scripts/checkm-log.R {params.log}  &> {log}
        Rscript panpiper/workflow/scripts/checkm_bin-stats.R {input.stats}  &> {log}
        """


rule fastani:
    input:
        ref=REFERENCE,
        file=os.path.join(OUT,"Quality/Assembly_filter/{file}/{file}.fna"),
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
        temp(os.path.join(OUT,"Quality/FastANI/{file}.txt")),
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
        stat=os.path.join(OUT,"Quality/CheckM/checkm_stats.txt"),
    params:
        checkm=os.path.join(OUT,"Quality/CheckM/checkm_log.txt"),
        outpath=os.path.join(OUT,"Quality"),
        ref=PARAMS_REF,
        ac=ANI_CUTOFF,
        cn=CONTIG_NUMBER,
        n=N50,
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
        python panpiper/workflow/scripts/filter_isolates.py -o {params.outpath} -a {input.ani} -l {params.checkm} -s {input.stat} -r {params.ref} -ac {params.ac} -cn {params.cn} -n {params.n} &> {log}
        """

rule print_results:
    input:
        ani=os.path.join(OUT,"Quality/FastANI/fastani_summary.txt"),
        stat=os.path.join(OUT,"Quality/CheckM/checkm_stats.txt"),
        log=os.path.join(OUT,"Quality/CheckM/checkm_log.txt"),
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
        mem="1G"
    threads: 1
    output:
        SAMPLES_OUT,
    shell:
        "Rscript -e \"rmarkdown::render('panpiper/workflow/scripts/quality_report.Rmd', output_dir = '{params.outdir}', params=list(checkm = '{input.stat}', log = '{input.log}', ani = '{input.ani}', passed = '{input.passed}', ref = '{params.ref}'))\" &> {log}"


