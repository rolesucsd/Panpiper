"""
Authors: Renee Oles
Date Updated: 1/4/2022
Purpose: Construct pangenome
"""
import os
import subprocess
import math
import distutils.util


OUT = config['out']
FASTA = config['fasta']
SAMPLE_LIST = config['list']
PARAMS = config['params']

with open(PARAMS, 'r') as fh:   
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

def read_filtered_files():
    with open(SAMPLE_LIST) as f:
        raw_reads = [sample for sample in f.read().split('\n') if len(sample) > 0]
        return(raw_reads)
filtered = read_filtered_files()

# Output
rule all:
    input:
        os.path.join(OUT,"Pangenome/AMR.txt"),
        os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.fna"),
        os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.bwt"),
        os.path.join(OUT,"Pangenome/Phylogroups/db/db_clusters.csv"),
        expand(os.path.join(OUT,"Pangenome/Prokka/{file}/{file}.gff"), file=filtered),
#        os.path.join(OUT,"Pangenome/Panaroo/eggnog.txt"),
#        expand(os.path.join(OUT,"Pangenome/Kraken/{file}.report"), file=filtered),


# Set up files for downstream analysis
rule setup:
    input:
        expand(os.path.join(FASTA, "{file}/{file}.fna"), file=filtered),
    params:
        os.path.join(OUT,"Pangenome")
    conda:
        "envs/amrfinder.yml"
    log:
        os.path.join(OUT,"report/setup.log"),
    output:
        os.path.join(OUT,"Pangenome/Unitig/unitig.list"),
        os.path.join(OUT,"Pangenome/Phylogroups/poppunk.list"),
    shell:
        """
        chmod u+x workflow/scripts/setup.sh
        workflow/scripts/setup.sh {input} {params} &> {log}
        amrfinder -u
        """


# Prokka will annotate individual assemblies and prepare them for pangenome construction
# Time requirement: between 5 and 10 minutes
rule prokka_multiple:
    input:
        gen=os.path.join(FASTA,"{file}/{file}.fna"),
    params:
        name="{file}",
        outdir=os.path.join(OUT,"Pangenome/Prokka/{file}"),
    conda:
        "envs/prokka.yml"
    log:
        os.path.join(OUT,"report/prokka_{file}.log"),
    output:
        gen=os.path.join(OUT,"Pangenome/Prokka/{file}/{file}.gff"),
        faa=os.path.join(OUT,"Pangenome/Prokka/{file}/{file}.faa"),
        fna=os.path.join(OUT,"Pangenome/Prokka/{file}/{file}.fna"),
    shell:
        "prokka --mincontiglen 500 --centre X --compliant --kingdom Bacteria --outdir {params.outdir} --prefix {params.name} --locustag {params.name} --force {input} --addgenes &> {log}"


# Panaroo is an updated version of Roary to create a pangenome
# 1 hour for 100 samples, scales linearlly
rule panaroo:
    input:
        gen=expand(os.path.join(OUT,"Pangenome/Prokka/{file}/{file}.gff"), file=filtered),
    params:
        os.path.join(OUT,"Pangenome/Panaroo"),
    log:
        os.path.join(OUT,"report/panaroo.log"),
    conda:
        "envs/panaroo.yml"
    output:
        os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.fa"),
        os.path.join(OUT,"Pangenome/Panaroo/core_gene_alignment.aln"),
    shell:
        "panaroo -i {input} -o {params} --remove-invalid-genes --clean-mode strict -a core --core_threshold 0.98 --len_dif_percent 0.98 -f 0.7 --merge_paralogs -t 20 &> {log}"


# Prokka will annotate the pangenome
rule prokka_pan:
    input:
        gen=os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.fa"),
    params:
        name="pan_genome_reference",
        outdir=os.path.join(OUT,"Pangenome/Panaroo/"),
    conda:
        "envs/prokka.yml"
    log:
        os.path.join(OUT,"report/prokka_pan.log"),
    output:
        gen=os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.gff"),
        faa=os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.faa"),
        fna=os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.fna"),
    shell:
        "prokka --mincontiglen 500 --centre X --compliant --kingdom Bacteria --outdir {params.outdir} --prefix {params.name} --locustag {params.name} --force {input} --addgenes &> {log}"


# The pangenome is indexed 
rule bwa_index_pan:
    input:
        os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.fa"),
    params:
        one=os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference"),
        inte=os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.fa.bwt"),
    conda:
        "envs/bwa.yml"
    log:
        os.path.join(OUT,"report/bwa_index_pan.log"),
    output:
        os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.bwt"),
    shell:
        """
        bwa index -a bwtsw {input} {params.one} &> {log}
        mv {params.inte} {output}
        """


# Create a starting tree based off core genome alignment
rule fasttree:
    input:
        os.path.join(OUT,"Pangenome/Panaroo/core_gene_alignment.aln"),
    conda:
        "envs/fasttree.yml"
    log:
        os.path.join(OUT,"report/fasttree.log"),
    output:
        os.path.join(OUT,"Pangenome/Phylogeny/fasttree.nwk"),
    shell:
        "FastTree -gtr -nt -gamma -seed 12345 {input} > {output} &> {log}"


# Build more accurate tree using fasttree as a starting point
rule raxml:
    input:
        tre=os.path.join(OUT,"Pangenome/Phylogeny/fasttree.nwk"),
        aln=os.path.join(OUT,"Pangenome/Panaroo/core_gene_alignment.aln"),
    params:
        outdir=os.path.join(OUT,"Pangenome/Phylogeny"),
    log:
        os.path.join(OUT,"report/raxml.log"),
    conda:
        "envs/raxml.yml"
    output:
        os.path.join(OUT,"Pangenome/Phylogeny/RAxML_result.tree2"),
    shell:
        "./standard-RAxML/raxmlHPC-AVX2 -m GTRCAT -F -f D -s {input.aln} -t {input.tre} -n tree2 -p 12345 -w {params.outdir} &> {log}"


# Continue to build tree from RAXML with iqtree
# all means combined tree search and bootstrapping analysis
# GTR indicates the sequence is DNA data
rule iqtree:
    input:
        aln=os.path.join(OUT,"Pangenome/Panaroo/core_gene_alignment.aln"),
        tre=os.path.join(OUT,"Pangenome/Phylogeny/RAxML_result.tree2"),
    conda:
        "envs/iqtree.yml"
    log:
        os.path.join(OUT,"report/iqtree.log"),
    output:
        os.path.join(OUT,"Pangenome/Panaroo/core_gene_alignment.aln.iqtree"),
    shell:
        "iqtree -m GTR+I+G -nt AUTO -s {input.aln} -te {input.tre} &> {log}"



# Identify antibiotic resistance in the all the samples
# TODO: Does snakemake allow optional parameters? 
rule amrfinder:
    input:
        setup=os.path.join(OUT,"Pangenome/Unitig/unitig.list"),
        fasta=os.path.join(OUT,"Pangenome/Prokka/{file}/{file}.faa"),
        gff=os.path.join(OUT,"Pangenome/Prokka/{file}/{file}.gff"),
    params:
        mut=os.path.join(OUT,"Pangenome/AMR/report_mut.txt"),
    conda:
        "envs/amrfinder.yml"
    log:
        os.path.join(OUT,"report/amrfinder_{file}.log"),
    output:
        os.path.join(OUT,"Pangenome/AMR/{file}.txt"),
    shell:
        """
        amrfinder -p {input.fasta} --plus --threads 20 --mutation_all {params.mut} -o {output} &> {log}
        """


# Summarize antibiotic resistance accross all files
rule concat_amr:
    input:
        expand(os.path.join(OUT,"Pangenome/AMR/{file}.txt"), file=filtered),
    log:
        os.path.join(OUT,"report/concat_amr.log"),
    output:
        os.path.join(OUT,"Pangenome/AMR.txt"),
    shell:
        "python scripts/concat_amrfinder.py {input} &> {log}"


# Mash- computes approximate distance between two samples QUICKLY
# locality-sensitive hashing (called MinHash sketches which have been used in lots of areas ie website comparisons)
# Time: 30 minutes for 571 samples
rule mash_fasta:
    input:
        expand(os.path.join(OUT,"Assembly/{file}/{file}.fna"), file=filtered),
    params:
        out=os.path.join(OUT,"Pangenome/Mash/fasta"),
        pref=os.path.join(OUT,"Assembly/Shovill/"),
    conda:
        "envs/mash.yml"
    log:
        os.path.join(OUT,"report/mash.log"),
    output:
        os.path.join(OUT,"Pangenome/Mash/fasta.tsv"),
    shell:
        """
        mash sketch -s 10000 -o {params.out} {input}
        mash dist {params.out}.msh {params.out}.msh -t > {output} &> {log}
        """

# Work on phylogroup division using Poppunk
rule poppunk_sketch:
    input:
        os.path.join(OUT,"Pangenome/Phylogroups/poppunk.list"),
    params:
        os.path.join(OUT,"Pangenome/Phylogroups/db"),
    conda:
        "envs/poppunk.yml"
    log:
        os.path.join(OUT,"report/poppunk_sketch.log"),
    output:
        os.path.join(OUT,"Pangenome/Phylogroups/db/db.dists.pkl"),
    shell:
        "poppunk --create-db --output {params} --r-files {input} --threads 20 &> {log}"


rule poppunk_qc:
    input:
        os.path.join(OUT,"Pangenome/Phylogroups/db/db.dists.pkl"),
    params:
        db=os.path.join(OUT,"Pangenome/Phylogroups/db"),
        ref=param_dict["ref"],
    conda:
        "envs/poppunk.yml"
    log:
        os.path.join(OUT,"report/poppunk_qc.log"),
    output:
        os.path.join(OUT,"Pangenome/Phylogroups/db/poppunk.txt"),
    shell:
        "poppunk --qc-db --ref-db {params.db} --type-isolate {params.ref} &> {log}"


rule poppunk_fit:
    input:
        os.path.join(OUT,"Pangenome/Phylogroups/db/db.dists.pkl"),
    params:
        db=os.path.join(OUT,"Pangenome/Phylogroups/db"),
        model="lineage",
    conda:
        "envs/poppunk.yml"
    log:
        os.path.join(OUT,"report/poppunk_fit.log"),
    output:
        os.path.join(OUT,"Pangenome/Phylogroups/db/db_clusters.csv"),
    shell:
        "poppunk --fit-model {params.model} --ref-db {params.db} --K 4 &> {log}"


rule virsorter_run:
    input:
        dbdir=directory("databases/virsorter"),
        fasta=os.path.join(FASTA,"{file}/{file}.fna"),
    params:
        outdir=os.path.join(OUT,"Pangenome/Virsorter/VS2_1/{file}"),
    conda:
        "envs/virsorter.yml"
    log:
        os.path.join(OUT,"report/virsorter_run1_{file}.log"),
    output:
        os.path.join(OUT,"Pangenome/Virsorter/VS2_1/{file}/final-viral-combined.fa"),
    shell:
        """g
        virsorter config --init-source --db-dir={params.db}
        virsorter run --keep-original-seq -i {input.fasta} -w {params.outdir} --include-groups dsDNAphage,ssDNA --min-length 5000 --min-score 0.5 -j 28 all &> {log}
        """


rule eggnog_mapper:
    input:
        fasta=os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.fa"),
    params:
        os.path.join(OUT,"Pangenome/Panaroo"),
    conda:
        "envs/eggnog.yml"
    log:
        os.path.join(OUT,"report/eggnog_mapper.log"),
    output:
        os.path.join(OUT,"Pangenome/Panaroo/eggnog.txt"),
    shell:
        """
        emapper.py -i {input.fasta} --data-dir databases/eggnog --sensmode more-sensitive -o {params} &> {log}
        """
