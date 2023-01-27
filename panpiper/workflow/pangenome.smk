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
        os.path.join(OUT,"Pangenome/Summary/amr.png"),
        os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.bwt"),
        os.path.join(OUT,"Pangenome/Summary/db_clusters.csv"),
        os.path.join(OUT,"Pangenome/Summary/core_gene_alignment.aln.iqtree"),
        os.path.join(OUT,"Pangenome/Summary/eggnog.txt"),
        os.path.join(OUT,"Pangenome/Summary/kraken_ag.txt"),
        os.path.join(OUT,"Pangenome/Unitig/unitig.pyseer"),
        os.path.join(OUT,"Pangenome/Summary/eggnog.txt"),


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
        chmod u+x panpiper/workflow/scripts/setup.sh
        panpiper/workflow/scripts/setup.sh {input} {params} &> {log}
        amrfinder -u
        """


# Bakta will annotate individual assemblies and prepare them for pangenome construction
# About 10 minutes per sample
rule bakta_multiple:
    input:
        gen=os.path.join(FASTA,"{file}/{file}.fna"),
    params:
        db=param_dict["bakta_dir"],
        name="{file}",
        outdir=os.path.join(OUT,"Pangenome/Bakta/{file}"),
    conda:
        "envs/bakta.yml"
    log:
        os.path.join(OUT,"report/bakta_{file}.log"),
    output:
        gen=os.path.join(OUT,"Pangenome/Bakta/{file}/{file}.gff3"),
        faa=os.path.join(OUT,"Pangenome/Bakta/{file}/{file}.faa"),
        fna=os.path.join(OUT,"Pangenome/Bakta/{file}/{file}.fna"),
    shell:
        "bakta --db {params.db} --output {params.outdir} --prefix {params.name} --locus-tag {params.name} {input} &> {log}"


# Panaroo is an updated version of Roary to create a pangenome
# 1 hour for 100 samples, scales linearlly
rule panaroo:
    input:
        gen=expand(os.path.join(OUT,"Pangenome/Bakta/{file}/{file}.gff3"), file=filtered),
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


# Bakta will annotate the pangenome
rule bakta_pan:
    input:
        gen=os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.fa"),
    params:
        db=param_dict["bakta_dir"],
        name="pan_genome_reference",
        outdir=os.path.join(OUT,"Pangenome/Summary/"),
    conda:
        "envs/bakta.yml"
    log:
        os.path.join(OUT,"report/bakta_pan.log"),
    output:
        gen=os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.gff3"),
        faa=os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.faa"),
        fna=os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.fna"),
    shell:
        "bakta --db {params.db} --output {params.outdir} --prefix {params.name} {input} &> {log}"

# Insert bakta pan (there is a bakta protein version that can annotate the pangenome)

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
        "FastTree -gtr -nt -gamma -seed 12345 {input} > {output}"


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
        "raxmlHPC-AVX2 -m GTRCAT -F -f D -s {input.aln} -t {input.tre} -p 12345 -n tree2 -w {params.outdir} &> {log}"


# Continue to build tree from RAXML with iqtree
# all means combined tree search and bootstrapping analysis
# GTR indicates the sequence is DNA data
rule iqtree:
    input:
        aln=os.path.join(OUT,"Pangenome/Panaroo/core_gene_alignment.aln"),
        tre=os.path.join(OUT,"Pangenome/Phylogeny/RAxML_result.tree2"),
    params:
        os.path.join(OUT,"Pangenome/Panaroo/core_gene_alignment.aln.iqtree"),
    conda:
        "envs/iqtree.yml"
    log:
        os.path.join(OUT,"report/iqtree.log"),
    output:
        os.path.join(OUT,"Pangenome/Summary/core_gene_alignment.aln.iqtree"),
    shell:
        """
        iqtree -m GTR+I+G -nt AUTO -s {input.aln} -te {input.tre} &> {log}
        cp {params.output} {output}
        """



# Identify antibiotic resistance in the all the samples
# TODO: Does snakemake allow optional parameters? 
rule amrfinder:
    input:
        setup=os.path.join(OUT,"Pangenome/Unitig/unitig.list"),
        fasta=os.path.join(OUT,"Pangenome/Bakta/{file}/{file}.faa"),
        gff=os.path.join(OUT,"Pangenome/Bakta/{file}/{file}.gff3"),
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
        os.path.join(OUT,"Pangenome/Summary/AMR.txt"),
    shell:
        "python panpiper/workflow/scripts/concat_amrfinder.py -i {input} -o {output} &> {log}"


# Summarize antibiotic resistance accross all files
rule graph_amr:
    input:
        os.path.join(OUT,"Pangenome/Summary/AMR.txt"),
    params:
        os.path.join(OUT,"Pangenome/Summary"),
    conda:
        "envs/r.yml"
    log:
        os.path.join(OUT,"report/graph_amr.log"),
    output:
        os.path.join(OUT,"Pangenome/Summary/amr.png"),
    shell:
        "Rscript panpiper/workflow/scripts/AMR_heatmap.R {input} {params} &> {log}"


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
        out=os.path.join(OUT,"Pangenome/Phylogroups/db/db_clusters.csv"),
        model="lineage",
    conda:
        "envs/poppunk.yml"
    log:
        os.path.join(OUT,"report/poppunk_fit.log"),
    output:
        os.path.join(OUT,"Pangenome/Summary/db_clusters.csv"),
    shell:
        """
        poppunk --fit-model {params.model} --ref-db {params.db} --K 4 &> {log}
        cp {params.out} {output}
        """


rule eggnog_mapper:
    input:
        protein=os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.faa"),
    params:
        outdir=os.path.join(OUT,"Pangenome/Summary"),
        db=param_dict['eggnog_dir'],
    conda:
        "envs/eggnog.yml"
    log:
        os.path.join(OUT,"report/eggnog_mapper.log"),
    output:
        os.path.join(OUT,"Pangenome/Summary/eggnog.txt"),
    shell:
        """
        emapper.py -i {input.protein} --data_dir {params.db} -o {params.outdir} --override &> {log}
        """

rule kraken:
    input:
        os.path.join(FASTA, "{file}/{file}.fna"),
    params:
        db=param_dict["kraken_dir"]
    conda:
        "envs/kraken.yml"
    log:
        os.path.join(OUT,"report/kraken_{file}.log"),
    output:
        out=os.path.join(OUT,"Pangenome/Kraken/{file}.out"),
        report=os.path.join(OUT,"Pangenome/Kraken/{file}.report")
    shell:
        """
        kraken2 --db {params.db} --output {output.out} --threads 20 --use-names --report {output.report} --use-mpa-style {input} &> {log}
        """

rule kraken_report:
    input:
        out=expand(os.path.join(OUT,"Pangenome/Kraken/{file}.out"), file=filtered),
    params:
        outpath=os.path.join(OUT,"Pangenome/Summary"),
    conda:
        "envs/r.yml"
    log:
        os.path.join(OUT,"report/kraken_summary.log"),
    output:
        report=os.path.join(OUT,"Pangenome/Summary/kraken_ag.txt")
    shell:
        """
        Rscript panpiper/workflow/scripts/kraken.R {input} {params.outpath} &> {log}
        """

# Unitig caller will create unitigs which is a kmer alternative (read about why using it in gwas is better or worse)
# In future note that you can use fastq or fasta or both see commented out line of code for other use
# Time: 25 minutes for 571 files 
rule unitig:
    input:
        os.path.join(OUT,"Pangenome/Unitig/unitig.list"),
    params:
        os.path.join(OUT,"Pangenome/Unitig/unitig"),
    conda:
        "envs/unitig.yml"
    log:
        os.path.join(OUT,"report/unitig.log"),
    output:
        os.path.join(OUT,"Pangenome/Unitig/unitig.pyseer"),
    shell:
        "unitig-caller --call --reads {input} --out {params} &> {log}" 
