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
Purpose: Construct pangenome
"""
import os
import sys
import subprocess
import math
import distutils.util
from scripts.wrangle_output import genes_long, gene_matrix_phylogroup, gene_matrix, process_data
from scripts.phylogroup_cluster import cluster_samples


OUT = config['out']
FASTA = config['fasta']
SAMPLE_LIST = config['list']
PARAMS = config['params']
PATH = config['scripts']

with open(PARAMS, 'r') as fh:   
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

BAKTA = param_dict["bakta_dir"]
REF = param_dict["ref"]
EGGNOG = param_dict['eggnog_dir']

def read_filtered_files():
    with open(SAMPLE_LIST) as f:
        raw_reads = [sample for sample in f.read().split('\n') if len(sample) > 0]
        return(raw_reads)
filtered = read_filtered_files()

# Output
rule all:
    input:
        os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.bwt"),
        os.path.join(OUT,"Pangenome/Summary/core_gene_alignment.aln.iqtree"),
        os.path.join(OUT,"Pangenome/Unitig/unitig.pyseer"),
        os.path.join(OUT,"Pangenome/Summary/genes_phylogroup_matrix_perc.txt"),
        os.path.join(OUT,"Pangenome/Summary/genes_anno.txt"),
        os.path.join(OUT,"Pangenome/Summary/genes_matrix.txt"),

# Download databases
rule bakta:
    input:
        BAKTA,
    conda:
        "envs/bakta.yml",
    log:
        os.path.join(OUT,"report/bakta_db.log"),
    benchmark:
        os.path.join(OUT,"benchmark/bakta.benchmark"),
    output:
        os.path.join(BAKTA, "db/bakta.db"),
    shell:
        "bakta_db download --output {input} &> {log}"

rule eggnog_db:
    input:
        EGGNOG,
    conda:
        "envs/eggnog.yml",
    log:
        os.path.join(OUT,"report/eggnog_db.log"),
    benchmark:
        os.path.join(OUT,"benchmark/eggnog_db.benchmark"),
    output:
        os.path.join(EGGNOG, "eggnog.db"),
    shell:
        "download_eggnog_data.py --data_dir {input}"

# Set up files for downstream analysis
rule setup:
    input:
        expand(os.path.join(FASTA, "{file}/{file}.fna"), file=filtered),
    params:
        outdir=os.path.join(OUT,"Pangenome"),
        path=PATH,
    conda:
        "envs/amrfinder.yml"
    log:
        os.path.join(OUT,"report/setup.log"),
    benchmark:
        os.path.join(OUT,"benchmark/setup.benchmark"),
    output:
        os.path.join(OUT,"Pangenome/Unitig/unitig.list"),
    shell:
        """
        chmod u+x {params.path}/setup.sh
        {params.path}/setup.sh {input} {params.outdir} &> {log}
        """

rule bakta_multiple:
    input:
        gen=os.path.join(FASTA,"{file}/{file}.fna"),
        db=os.path.join(BAKTA, "db/bakta.db"),
    params:
        db=os.path.join(BAKTA, "db"),
        name="{file}",
        outdir=os.path.join(OUT,"Pangenome/Bakta/{file}"),
    conda:
        "envs/bakta.yml"
    log:
        os.path.join(OUT,"report/bakta_{file}.log"),
    benchmark:
        os.path.join(OUT,"benchmark/bakta_{file}.benchmark"),
    output:
        gen=os.path.join(OUT,"Pangenome/Bakta/{file}/{file}.gff3"),
        faa=os.path.join(OUT,"Pangenome/Bakta/{file}/{file}.faa"),
        fna=os.path.join(OUT,"Pangenome/Bakta/{file}/{file}.fna"),
    shell:
#        "bakta --skip-plot --db {params.db} --output {params.outdir} --prefix {params.name} --locus-tag {params.name} {input.gen} &> {log}"
        """
        bakta --skip-plot --db {params.db} --output {params.outdir} --prefix {params.name} {input.gen} &> {log}
        """


# Panaroo is an updated version of Roary to create a pangenome
# 1 hour for 100 samples, scales linearlly
rule panaroo:
    input:
        gen=expand(os.path.join(OUT,"Pangenome/Bakta/{file}/{file}.gff3"), file=filtered),
    params:
        os.path.join(OUT,"Pangenome/Panaroo"),
    conda:
        "envs/panaroo.yml"
    log:
        os.path.join(OUT,"report/panaroo.log"),
    benchmark:
        os.path.join(OUT,"benchmark/panaroo.benchmark"),
    output:
        os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.fa"),
        os.path.join(OUT,"Pangenome/Panaroo/core_gene_alignment.aln"),
        os.path.join(OUT,"Pangenome/Panaroo/gene_presence_absence_roary.csv"),
    shell:
        "panaroo -i {input} -o {params} --remove-invalid-genes --clean-mode strict -a core --core_threshold 0.98 --len_dif_percent 0.98 -f 0.7 --merge_paralogs -t 20 &> {log}"

# Translate DNA to AA
rule translate:
    input:
        os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.fa"),
    params:
        path=PATH,
    output:
        gen=os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.faa"),
    shell:
        "python {params.path}/nc_aa_translate.py {input} {output}"

# Bakta will annotate the pangenome
rule bakta_pan:
    input:
        gen=os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.fa"),
        db=os.path.join(BAKTA, "db/pfam.h3p"),
    params:
        db=os.path.join(BAKTA, "db"),
        name="pan_genome_reference",
        outdir=os.path.join(OUT,"Pangenome/Summary/"),
    conda:
        "envs/bakta.yml"
    log:
        os.path.join(OUT,"report/bakta_pan.log"),
    benchmark:
        os.path.join(OUT,"benchmark/bakta_pan.benchmark"),
    output:
        genes=os.path.join(OUT,"Pangenome/Summary/pan_genome_reference.tsv"),
        fasta=os.path.join(OUT,"Pangenome/Summary/pan_genome_reference.faa"),
        gff=os.path.join(OUT,"Pangenome/Summary/pan_genome_reference.gff3"),
    shell:
        "bakta --skip-plot --db {params.db} --output {params.outdir} --prefix {params.name} {input.gen} &> {log}"

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
    benchmark:
        os.path.join(OUT,"benchmark/bwa_index_pan.benchmark"),
    output:
        os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.bwt"),
    shell:
        """
        bwa index -a bwtsw {input} {params.one} &> {log}
        mv {params.inte} {output}
        """

# Create a starting tree based off core genome alignment
# gtr nt = a tree for a nucleotide alignment with GTR+CAT model 
# gamma = 5% slower but rescales teh branch lengthsa nd computes a Gamma20-based likelihood which is more comparable across runs 
rule fasttree:
    input:
        os.path.join(OUT,"Pangenome/Panaroo/core_gene_alignment.aln"),
    conda:
        "envs/fasttree.yml"
    log:
        os.path.join(OUT,"report/fasttree.log"),
    benchmark:
        os.path.join(OUT,"benchmark/fasttree.benchmark"),
    output:
        os.path.join(OUT,"Pangenome/Phylogeny/fasttree.nwk"),
    shell:
        "FastTree  -intree -gtr -nt -gamma -seed 12345 {input} > {output}"


# Build more accurate tree using fasttree as a starting point
# This step infers topology 
# See what these parameters meant in old method: and GTRCAT -F -f D
rule raxml:
    input:
        tre=os.path.join(OUT,"Pangenome/Phylogeny/fasttree.nwk"),
        aln=os.path.join(OUT,"Pangenome/Panaroo/core_gene_alignment.aln"),
    params:
        outdir=os.path.join(OUT,"Pangenome/Phylogeny/two"),
    conda:
        "envs/raxml.yml",
    log:
        os.path.join(OUT,"report/raxml.log"),
    benchmark:
        os.path.join(OUT,"benchmark/raxml.benchmark"),
    output:
        os.path.join(OUT,"Pangenome/Phylogeny/two.raxml.bestTree"),
    shell:
        "raxml-ng --model GTR+G --msa {input.aln} --tree {input.tre} --seed 12345 --prefix {params.outdir} &> {log}"


# Continue to build tree from RAXML with iqtree
# all means combined tree search and bootstrapping analysis
# GTR indicates the sequence is DNA data
# This step optimizes branch lengths and computes likelihood score 
rule iqtree:
    input:
        aln=os.path.join(OUT,"Pangenome/Panaroo/core_gene_alignment.aln"),
        tre=os.path.join(OUT,"Pangenome/Phylogeny/two.raxml.bestTree"),
    params:
        output=os.path.join(OUT,"Pangenome/Panaroo/core_gene_alignment.aln.iqtree"),
    conda:
        "envs/iqtree.yml"
    log:
        os.path.join(OUT,"report/iqtree.log"),
    benchmark:
        os.path.join(OUT,"benchmark/iqtree.benchmark"),
    output:
        os.path.join(OUT,"Pangenome/Summary/core_gene_alignment.aln.iqtree"),
    shell:
        """
        iqtree -m GTR+I+G -nt AUTO -s {input.aln} -te {input.tre} &> {log}
        cp {params.output} {output}
        """

rule eggnog_mapper:
    input:
        protein=os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.faa"),
        db=os.path.join(EGGNOG, "eggnog.db"),
        fna=os.path.join(OUT,"Pangenome/Summary/pan_genome_reference.tsv"),
    params:
        outdir=os.path.join(OUT,"Pangenome/Summary/Summary"),
        db=EGGNOG,
    conda:
        "envs/eggnog.yml"
    log:
        os.path.join(OUT,"report/eggnog_mapper.log"),
    benchmark:
        os.path.join(OUT,"benchmark/eggnog_mapper.benchmark"),
    output:
        os.path.join(OUT,"Pangenome/Summary/Summary.emapper.annotations"),
    shell:
        """
        emapper.py -i {input.protein} --data_dir {params.db} -o {params.outdir} --override &> {log}
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
    benchmark:
        os.path.join(OUT,"benchmark/unitig.benchmark"),
    output:
        os.path.join(OUT,"Pangenome/Unitig/unitig.pyseer"),
    shell:
        "unitig-caller --call --reads {input} --out {params} &> {log}" 

rule mash_fasta:
    input:
        expand(os.path.join(FASTA, "{file}/{file}.fna"), file=filtered),
    params:
        os.path.join(OUT,"Pangenome/Mash/fasta"),
    conda:
        "envs/mash.yml"
    output:
        os.path.join(OUT,"Pangenome/Mash/fasta.tsv"),
    shell:
        """
        mash sketch -s 10000 -o {params} {input}
        mash dist {params}.msh {params}.msh -t > {output}
        """

# Add rule to interpret mash output
rule phylogroups:
    input:
        os.path.join(OUT,"Pangenome/Mash/fasta.tsv"),
    params:
        os.path.join(OUT,"Pangenome/Mash/"),
    output:
        os.path.join(OUT,"Pangenome/Mash/phylogroups.txt"),
    run:
        cluster_samples(mash_file=os.path.join(OUT,"Pangenome/Mash/fasta.tsv"), output_folder=os.path.join(OUT,"Pangenome/Mash/"))


# Add rule to wrangle panaroo output 
rule process_data:
    input:
        eggnog=os.path.join(OUT,"Pangenome/Summary/Summary.emapper.annotations"),
        bakta=os.path.join(OUT,"Pangenome/Summary/pan_genome_reference.tsv"),
    params:
        output=os.path.join(OUT,"Pangenome/Summary/"),
        genes=os.path.join(OUT,"Pangenome/Panaroo/gene_presence_absence_roary.csv"),
    output:
        os.path.join(OUT,"Pangenome/Summary/genes_anno.txt"),
    run:
        process_data(bakta=os.path.join(OUT,"Pangenome/Summary/pan_genome_reference.tsv"), genes=os.path.join(OUT,"Pangenome/Panaroo/gene_presence_absence_roary.csv"), eggnog=os.path.join(OUT,"Pangenome/Summary/Summary.emapper.annotations"), output=os.path.join(OUT,"Pangenome/Summary/"))

rule genes_long:
    input:
        clusters=os.path.join(OUT,"Pangenome/Mash/phylogroups.txt"),
        genes=os.path.join(OUT,"Pangenome/Panaroo/gene_presence_absence_roary.csv"),
    params:
        output=os.path.join(OUT,"Pangenome/Summary/"),
    output:
        os.path.join(OUT,"Pangenome/Summary/genes_long.txt"),
    run:
        genes_long(clusters=os.path.join(OUT,"Pangenome/Mash/phylogroups.txt"), genes=os.path.join(OUT,"Pangenome/Panaroo/gene_presence_absence_roary.csv"), output=os.path.join(OUT,"Pangenome/Summary/"))


rule gene_matrix_phylogroup:
    input:
        genes_long=os.path.join(OUT,"Pangenome/Summary/genes_long.txt"),
    params:
        output=os.path.join(OUT,"Pangenome/Summary/"),
    output:
        os.path.join(OUT,"Pangenome/Summary/genes_phylogroup_matrix_perc.txt"),
    run:
        gene_matrix_phylogroup(genes_long=os.path.join(OUT,"Pangenome/Summary/genes_long.txt"), output=os.path.join(OUT,"Pangenome/Summary/"))

rule gene_matrix:
    input:
        genes=os.path.join(OUT,"Pangenome/Panaroo/gene_presence_absence_roary.csv"),
    params:
        output=os.path.join(OUT,"Pangenome/Summary/"),
    output:
        os.path.join(OUT,"Pangenome/Summary/genes_matrix.txt"),
    run:
        gene_matrix(genes=os.path.join(OUT,"Pangenome/Panaroo/gene_presence_absence_roary.csv"), output=os.path.join(OUT,"Pangenome/Summary/"))

