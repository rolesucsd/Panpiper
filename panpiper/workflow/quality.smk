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
SAMPLES_OUT = os.path.join(OUT, 'Quality/"quality_report.html"')

with open(PARAMS, 'r') as fh:   
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

GC = param_dict["gc"]
GENOME_SIZE = param_dict["genome_size"]
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
        assembly=os.path.join(FASTA,"{file}/contigs.fa"),
    params:
        contig=500,
    conda:
        "envs/bbmap.yml"
    log:
        os.path.join(OUT,"report/prokka_filter_{file}.log"),
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
        threads=40,
        binset=os.path.join(OUT,"Quality/Assembly_filter/{file}"),
    output:
        stats=os.path.join(OUT,"Quality/Assembly_filter/{file}/storage/bin_stats.analyze.tsv"),
    conda:
        "envs/checkm.yml"
    log:
        log=os.path.join(OUT,"Quality/Assembly_filter/{file}/lineage.log"),
    shell:
#        echo ".snakemake/conda/checkm" | checkm data setRoot "../.snakemake/conda/checkm"
        """
        checkm lineage_wf -t 20 -x fna {params.binset} {params.binset} &> {log}
        """


rule checkm_to_graph:
    input:
        stats=expand(os.path.join(OUT,"Quality/Assembly_filter/{file}/storage/bin_stats.analyze.tsv"), file=READS),
    params:
        log=expand(os.path.join(OUT,"Quality/Assembly_filter/{file}/lineage.log"), file=READS),
    output:
        png=os.path.join(OUT,"Quality/CheckM/checkm_log.txt"),
        stats=os.path.join(OUT,"Quality/CheckM/checkm_stats.txt"),
    conda:
        "envs/r.yml"
    log:
        os.path.join(OUT,"report/checkm_graph.log"),
    shell:
        """
        Rscript panpiper/workflow/scripts/checkm-log.R {params.log}  &> {log}
        Rscript panpiper/workflow/scripts/checkm_bin-stats.R {input.stats}  &> {log}
        """


rule fastani_list_create:
    input:
        fasta=expand(os.path.join(OUT,"Quality/Assembly_filter/{file}/{file}.fna"), file=READS),
    params:
        ref=REFERENCE,
    output:
        os.path.join(OUT,"Quality/FastANI/fastani_list.txt"),
    log:
        os.path.join(OUT,"report/fastani_list_create.log"),
    shell:
        """
        chmod u+x panpiper/workflow/scripts/create_list.sh
        panpiper/workflow/scripts/create_list.sh {input} {params} {output} &> {log}
        """


rule fastani:
    input:
        os.path.join(OUT,"Quality/FastANI/fastani_list.txt"),
    output:
        os.path.join(OUT,"Quality/FastANI/matrix.txt"),
    conda:
        "envs/fastani.yml"
    log:
        os.path.join(OUT,"report/fastani.log"),
    shell:
        "fastANI --ql {input} --rl {input} -o {output} -t 20 &> {log}"


rule fastani_long_to_wide:
    input:
        os.path.join(OUT,"Quality/FastANI/matrix.txt"),
    output:
        txt=os.path.join(OUT,"Quality/FastANI/matrix_wide.txt"),
    log:
        os.path.join(OUT,"report/fastani_long_to_wide.log"),
    shell:
        "python panpiper/workflow/scripts/long_to_wide.py {input} {output.txt} &> {log}"


# Filter files based off user-defined rules
# Need to edit this python file bc it's messy rn
rule filter_files:
    input:
        ani=os.path.join(OUT,"Quality/FastANI/matrix_wide.txt"),
        stat=os.path.join(OUT,"Quality/CheckM/checkm_stats.txt"),
    params:
        checkm=os.path.join(OUT,"Quality/CheckM/checkm_log.txt"),
        outpath=os.path.join(OUT,"Quality"),
        ref=PARAMS_REF,
        gc=GC,
        genome_size=GENOME_SIZE,
        ac=ANI_CUTOFF,
        cn=CONTIG_NUMBER,
        n=N50,
    log:
        os.path.join(OUT,"report/filter_files.log"),
    output:
        os.path.join(OUT,"Quality/sample_list.txt"),
    shell:
        """
        python panpiper/workflow/scripts/filter_isolates.py -o {params.outpath} -a {input.ani} -l {params.checkm} -s {input.stat} -r {params.ref} -gc {params.gc} -g {params.genome_size} -ac {params.ac} -cn {params.cn} -n {params.n} &> {log}
        """

rule print_results:
    input:
        ani=os.path.join(OUT,"Quality/FastANI/matrix_wide.txt"),
        stat=os.path.join(OUT,"Quality/CheckM/checkm_stats.txt"),
        log=os.path.join(OUT,"Quality/CheckM/checkm_log.txt"),
        passed=os.path.join(OUT,"Quality/sample_list.txt"),
    params:
        ref=PARAMS_REF,
    conda:
        "envs/r.yml"
    log:
        os.path.join(OUT,"report/print_files.log"),
    output:
        SAMPLES_OUT,
    shell:
        "Rscript -e \"rmarkdown::render('panpiper/workflow/scripts/quality_report.Rmd', params=list(checkm = '{input.stat}', log = '{input.log}', ani = '{input.ani}', passed = '{input.passed}', ref = '{params.ref}'))\" &> {log}"

