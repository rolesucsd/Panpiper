"""
Authors: Renee Oles
Date Updated: 8/5/21
Purpose: Snakemake file that takes assemblies and performs quality control and filtering
Input: 
    Isolate reference(s)- fasta file
    Sample assemblies - fasta files
Example usage: snakemake -s assembly -j 10 --use-conda --rerun-incomplete --cluster-config cluster.json --cluster "sbatch -A {cluster.account} --mem {cluster.mem} -t {cluster.time} --cpus-per-task {cluster.cpus}" &
Ouput: CheckM files, ANI, and txt file of samples passing all filters 
"""

(reads,) = glob_wildcards("resources/fasta/{file}/contigs.fa")


# Output
rule all:
    input:
        "results/Quality/sample_list.txt",


# Remove contigs smaller than 500 bp
# 30 sec per file
rule contig_filter:
    input:
        assembly="resources/fasta/{file}/contigs.fa",
    params:
        contig=500,
    conda:
        "../envs/quality.yml"
    log:
        "workflow/report/prokka_filter_{file}.log",
    output:
        fna="results/Assembly/{file}/{file}.fna",
    shell:
        "reformat.sh in={input.assembly} out={output} minlength={params.contig} &> {log}"


# 5 minutes per sample
rule run_checkm:
    input:
        file="results/Assembly/{file}/{file}.fna",
        binset="results/Assembly/{file}",
    params:
        threads=40,
        prefix="results/Quality/Checkm/{file}",
    output:
        stats="results/Quality/Checkm/{file}/storage/bin_stats.analyze.tsv",
    conda:
        "../envs/quality.yml"
    log:
        log="results/Quality/Checkm/{file}/lineage.log",
    shell:
        "checkm lineage_wf -t 20 -x fna {input.binset} {params.prefix} &> {log}"


rule checkm_to_graph:
    input:
        stats=expand(
            "results/Quality/Checkm/{file}/storage/bin_stats.analyze.tsv", file=reads
        ),
    params:
        log=expand("results/Quality/Checkm/{file}/lineage.log", file=reads),
    output:
        png="results/Quality/checkm_log.txt",
        stats="results/Quality/checkm_stats.txt",
    conda:
        "../envs/quality.yml"
    log:
        "workflow/report/checkm_graph.log",
    shell:
        """
        Rscript workflow/scripts/checkm-log.R {params.log}  &> {log}
        Rscript workflow/scripts/checkm_bin-stats.R {input.stats}  &> {log}
        """


rule fastani_list_create:
    input:
        fasta=expand("results/Assembly/{file}/{file}.fna", file=reads),
    params:
        ref=config["reference_fna"],
    output:
        "results/Quality/fastani_list.txt",
    log:
        "workflow/report/fastani_list_create.log",
    shell:
        """
        chmod u+x workflow/scripts/create_list.sh
        workflow/scripts/create_list.sh {input} {params} &> {log}
        """


rule fastani:
    input:
        "results/Quality/fastani_list.txt",
    output:
        "results/Quality/matrix.txt",
    conda:
        "../envs/quality.yml"
    log:
        "workflow/report/fastani.log",
    shell:
        "fastANI --ql {input} --rl {input} -o {output} -t 20 &> {log}"


rule fastani_long_to_wide:
    input:
        "results/Quality/matrix.txt",
    output:
        txt="results/Quality/matrix_wide.txt",
    log:
        "workflow/report/fastani_long_to_wide.log",
    shell:
        "python workflow/scripts/long_to_wide.py {input} {output.txt} &> {log}"


# Filter files based off user-defined rules
rule filter_files:
    input:
        ani="results/Quality/matrix_wide.txt",
        stat="results/Quality/checkm_stats.txt",
    params:
        checkm="results/Quality/checkm_log.txt",
        ref=config["ref"],
        gc=config["gc"],
        genome_size=config["genome_size"],
    log:
        "workflow/report/filter_files.log",
    output:
        "results/Quality/sample_list.txt",
    shell:
        """
        python workflow/scripts/filter_isolates.py -o {output} -a {input.ani} -l {params.checkm} -s {input.stat} -r {params.ref} -gc {params.gc} -g {params.genome_size} &> {log}
        chmod u+x workflow/scripts/copy_files.sh
        workflow/scripts/copy_files.sh
        """