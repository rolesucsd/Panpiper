"""
Authors: Renee Oles
Date Updated: 11/1/22
Purpose: Mgwas snakemake
Config file: 
    fasta_path - path to fasta files
    fastq_path - path to fastq files
    reference - path to reference file in gbk format
Input: 
Example usage: nohup snakemake -s mgwas.smk -j 20 --use-conda --cluster-config cluster.json --cluster "sbatch -A {cluster.account} --mem {cluster.mem} -t {cluster.time} --cpus-per-task {cluster.cpus}" &
Output:         
"""

# find . -type f -name "*_genomic.fna" -printf "/%P\n" | while read FILE ; do DIR=$(dirname "$FILE" );\
# mv ."$FILE" ."$DIR""$DIR".fna;
# done

(filtered,) = glob_wildcards("resources/fasta/{file}.fna")


rule all:
    input:
        "results/Pangenome/AMR.txt",
        "results/Pangenome/Panaroo/pan_genome_reference.fna",
        "results/Pangenome/Panaroo/pan_genome_reference.bwt",
        "results/Pangenome/Phylogroups/db/db_clusters.csv",
        expand("results/Pangenome/Prokka/{file}/{file}.gff", file=filtered),


rule setup:
    input:
        "results/Quality/sample_list.txt",
    output:
        "results/Pangenome/Unitig/unitig.list",
        "results/Pangenome/Fsm-lite/fsm-lite.list",
        "results/Pangenome/Phylogroups/poppunk.list",
    shell:
        """
        chmod u+x scripts/setup.sh
        scripts/setup.sh {input}
        """


# Time requirement: between 5 and 10 minutes
rule prokka_multiple:
    input:
        gen="results/Assembly/{file}/{file}.fna",
    params:
        name="{file}",
        outdir="results/Pangenome/Prokka/{file}",
    conda:
        "../envs/pangenome.yml"
    log:
        "log/prokka_{file}.log",
    output:
        gen="results/Pangenome/Prokka/{file}/{file}.gff",
        faa="results/Pangenome/Prokka/{file}/{file}.faa",
        fna="results/Pangenome/Prokka/{file}/{file}.fna",
    shell:
        "prokka --mincontiglen 500 --centre X --compliant --kingdom Bacteria --outdir {params.outdir} --prefix {params.name} --locustag {params.name} --force {input} --addgenes &> {log}"


# 1 hour for 100 samples, scales linearlly
rule panaroo:
    input:
        gen=expand("results/Pangenome/Prokka/{file}/{file}.gff", file=filtered),
    params:
        "results/Pangenome/Panaroo",
    log:
        "log/panaroo.log",
    output:
        "results/Pangenome/Panaroo/pan_genome_reference.fa",
        "results/Pangenome/Panaroo/core_gene_alignment.aln",
    shell:
        "panaroo -i {input} -o {params} --remove-invalid-genes --clean-mode strict -a core --core_threshold 0.98 --len_dif_percent 0.98 -f 0.7 --merge_paralogs -t 20 &> {log}"


rule prokka_pan:
    input:
        gen="results/Pangenome/Panaroo/pan_genome_reference.fa",
    params:
        name="pan_genome_reference",
        outdir="results/Pangenome/Panaroo/",
    conda:
        "../envs/pangenome.yml"
    log:
        "log/prokka_pan.log",
    output:
        gen="results/Pangenome/Panaroo/pan_genome_reference.gff",
        faa="results/Pangenome/Panaroo/pan_genome_reference.faa",
        fna="results/Pangenome/Panaroo/pan_genome_reference.fna",
    shell:
        "prokka --mincontiglen 500 --centre X --compliant --kingdom Bacteria --outdir {params.outdir} --prefix {params.name} --locustag {params.name} --force {input} --addgenes &> {log}"


rule bwa_index_pan:
    input:
        "results/Pangenome/Panaroo/pan_genome_reference.fa",
    params:
        one="results/Pangenome/Panaroo/pan_genome_reference",
        inte="results/Pangenome/Panaroo/pan_genome_reference.fa.bwt",
    conda:
        "../envs/pangenome.yml"
    output:
        "results/Pangenome/Panaroo/pan_genome_reference.bwt",
    shell:
        """
        bwa index -a bwtsw {input} {params.one}
        mv {params.inte} {output}
        """


rule fasttree:
    input:
        "results/Pangenome/Panaroo/core_gene_alignment.aln",
    conda:
        "../envs/pangenome.yml"
    output:
        "results/Pangenome/Phylogeny/fasttree.nwk",
    shell:
        "FastTree -gtr -nt -gamma -seed 12345 {input} > {output}"


rule raxml:
    input:
        tre="results/Pangenome/Phylogeny/fasttree.nwk",
        aln="results/Pangenome/Panaroo/core_gene_alignment.aln",
    params:
        outdir="results/Pangenome/Phylogeny",
    output:
        "results/Pangenome/Phylogeny/RAxML_result.tree2",
    shell:
        "./standard-RAxML/raxmlHPC-AVX2 -m GTRCAT -F -f D -s {input.aln} -t {input.tre} -n tree2 -p 12345 -w {params.outdir}"


# all means combined tree search and bootstrapping analysis
# GTR indicates the sequence is DNA data


rule iqtree:
    input:
        aln="results/Pangenome/Panaroo/core_gene_alignment.aln",
        tre="results/Pangenome/Phylogeny/RAxML_result.tree2",
    conda:
        "../envs/pangenome.yml"
    log:
        "log/iqtree.log",
    output:
        "results/Pangenome/Panaroo/core_gene_alignment.aln.iqtree",
    shell:
        "iqtree -m GTR+I+G -nt AUTO -s {input.aln} -te {input.tre}"


rule amrfinder:
    input:
        fasta="results/Pangenome/Prokka/{file}/{file}.faa",
        gff="results/Pangenome/Prokka/{file}/{file}.gff",
    params:
        mut="results/Pangenome/AMR/report_mut.txt",
    conda:
        "../envs/pangenome.yml"
    output:
        "results/Pangenome/AMR/{file}.txt",
    shell:
        "amrfinder -d /home/roles/anaconda3/../envs/amrfinder/share/amrfinderplus/data/2022-05-26.1 -p {input.fasta} --plus --threads 20 --mutation_all {params.mut} -o {output}"


rule concat_amr:
    input:
        expand("results/Pangenome/AMR/{file}.txt", file=filtered),
    output:
        "results/Pangenome/AMR.txt",
    shell:
        "python scripts/concat_amrfinder.py {input}"


# Mash- computes approximate distance between two samples QUICKLY
# locality-sensitive hashing (called MinHash sketches which have been used in lots of areas ie website comparisons)
# Time: 30 minutes for 571 samples
rule mash_fasta:
    input:
        expand("results/Assembly/{file}/{file}.fna", file=filtered),
    params:
        out="results/Pangenome/Mash/fasta",
        pref="results/Assembly/Shovill/",
    conda:
        "../envs/pangenome.yml"
    output:
        "results/Pangenome/Mash/fasta.tsv",
    shell:
        """
        mash sketch -s 10000 -o {params.out} {input}
        mash dist {params.out}.msh {params.out}.msh -t > {output}
        """

rule poppunk_sketch:
    input:
        "results/Pangenome/Phylogroups/poppunk.list",
    params:
        "results/Pangenome/Phylogroups/db",
    conda:
        "../envs/pangenome.yml"
    output:
        "results/Pangenome/Phylogroups/db/db.dists.pkl",
    shell:
        "poppunk --create-db --output {params} --r-files {input} --threads 20"


rule poppunk_qc:
    input:
        "results/Pangenome/Phylogroups/db/db.dists.pkl",
    params:
        db="results/Pangenome/Phylogroups/db",
        ref=config["ref"],
    conda:
        "../envs/pangenome.yml"
    output:
        "results/Pangenome/Phylogroups/db/popunk.txt",
    shell:
        "poppunk --qc-db --ref-db {params.db} --type-isolate {params.ref}"


rule poppunk_fit:
    input:
        "results/Pangenome/Phylogroups/db/db.dists.pkl",
    params:
        db="results/Pangenome/Phylogroups/db",
        model="lineage",
    conda:
        "../envs/pangenome.yml"
    output:
        "results/Pangenome/Phylogroups/db/db_clusters.csv",
    shell:
        "poppunk --fit-model {params.model} --ref-db {params.db} --K 4"