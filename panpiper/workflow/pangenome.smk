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

BAKTA = param_dict["bakta_dir"]
REF = param_dict["ref"]
EGGNOG = param_dict['eggnog_dir']
KRAKEN = param_dict["kraken_dir"]

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
        os.path.join(OUT,"Pangenome/Summary/Summary.emapper.annotations"),
#        os.path.join(OUT,"Pangenome/Summary/kraken_ag.txt"),
        os.path.join(OUT,"Pangenome/Unitig/unitig.pyseer"),

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
        os.path.join(BAKTA, "bakta.db"),
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

rule kraken_db:
    input:
        KRAKEN,
    conda:
        "envs/kraken.yml",
    log:
        os.path.join(OUT,"report/kraken_db.log"),
    benchmark:
        os.path.join(OUT,"benchmark/kraken_db.benchmark"),
    output:
        os.path.join(KRAKEN, "taxo.k2d"),
    shell:
        "kraken2-build --standard -db {input}"

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
    benchmark:
        os.path.join(OUT,"benchmark/setup.benchmark"),
    output:
        os.path.join(OUT,"Pangenome/Unitig/unitig.list"),
        os.path.join(OUT,"Pangenome/Phylogroups/poppunk.list"),
    shell:
        """
        chmod u+x panpiper/workflow/scripts/setup.sh
        panpiper/workflow/scripts/setup.sh {input} {params} &> {log}
        """


# Bakta will annotate individual assemblies and prepare them for pangenome construction
# About 10 minutes per sample
rule bakta_multiple:
    input:
        gen=os.path.join(FASTA,"{file}/{file}.fna"),
        db=os.path.join(BAKTA, "bakta.db"),
    params:
        db=BAKTA,
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
        "bakta --skip-plot --db {params.db} --output {params.outdir} --prefix {params.name} {input.gen} &> {log}"


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
    shell:
        "panaroo -i {input} -o {params} --remove-invalid-genes --clean-mode strict -a core --core_threshold 0.98 --len_dif_percent 0.98 -f 0.7 --merge_paralogs -t 20 &> {log}"

# Translate DNA to AA
rule translate:
    input:
        os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.fa"),
    output:
        gen=os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.faa"),
    shell:
        "python panpiper/workflow/scripts/nc_aa_translate.py {input} {output}"

# Bakta will annotate the pangenome
rule bakta_pan:
    input:
        gen=os.path.join(OUT,"Pangenome/Panaroo/pan_genome_reference.faa"),
        db=os.path.join(BAKTA, "bakta.db"),
    params:
        db=BAKTA,
        name="pan_genome_reference",
        outdir=os.path.join(OUT,"Pangenome/Summary/"),
    conda:
        "envs/bakta.yml"
    log:
        os.path.join(OUT,"report/bakta_pan.log"),
    benchmark:
        os.path.join(OUT,"benchmark/bakta_pan.benchmark"),
    output:
        fna=os.path.join(OUT,"Pangenome/Summary/pan_genome_reference.tsv"),
    shell:
        "bakta_proteins --db {params.db} --output {params.outdir} --prefix {params.name} {input.gen} &> {log}"

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
        "FastTree -gtr -nt -gamma -seed 12345 {input} > {output}"


# Build more accurate tree using fasttree as a starting point
rule raxml:
    input:
        tre=os.path.join(OUT,"Pangenome/Phylogeny/fasttree.nwk"),
        aln=os.path.join(OUT,"Pangenome/Panaroo/core_gene_alignment.aln"),
    params:
        outdir=os.path.join(OUT,"Pangenome/Phylogeny"),
    conda:
        "envs/raxml.yml",
    log:
        os.path.join(OUT,"report/raxml.log"),
    benchmark:
        os.path.join(OUT,"benchmark/raxml.benchmark"),
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

# Identify antibiotic resistance in the all the samples
# TODO: Does snakemake allow optional parameters? 
rule amrfinder:
    input:
        fasta=os.path.join(OUT,"Pangenome/Bakta/{file}/{file}.faa"),
        gff=os.path.join(OUT,"Pangenome/Bakta/{file}/{file}.gff3"),
    params:
        mut=os.path.join(OUT,"Pangenome/AMR/report_mut.txt"),
    conda:
        "envs/amrfinder.yml"
    log:
        os.path.join(OUT,"report/amrfinder_{file}.log"),
    benchmark:
        os.path.join(OUT,"benchmark/armfinder_{file}.benchmark"),
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
    resources:
        mem="1G"
    threads: 1
    benchmark:
        os.path.join(OUT,"benchmark/concat_amr.benchmark"),
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
    benchmark:
        os.path.join(OUT,"benchmark/graph_amr.benchmark"),
    resources:
        mem="1G"
    threads: 1
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
    benchmark:
        os.path.join(OUT,"benchmark/poppunk_sketch.benchmark"),
    resources:
        mem="5G"
    output:
        os.path.join(OUT,"Pangenome/Phylogroups/db/db.dists.pkl"),
    shell:
        "poppunk --create-db --output {params} --r-files {input} --threads 20 &> {log}"


rule poppunk_qc:
    input:
        os.path.join(OUT,"Pangenome/Phylogroups/db/db.dists.pkl"),
    params:
        db=os.path.join(OUT,"Pangenome/Phylogroups/db"),
        ref=REF,
    conda:
        "envs/poppunk.yml"
    log:
        os.path.join(OUT,"report/poppunk_qc.log"),
    benchmark:
        os.path.join(OUT,"benchmark/poppunk_qc.benchmark"),
    resources:
        mem="5G"
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
    benchmark:
        os.path.join(OUT,"benchmark/poppunk_fit.benchmark"),
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

rule kraken:
    input:
        file=os.path.join(FASTA, "{file}/{file}.fna"),
        db=os.path.join(KRAKEN, "taxo.k2d"),
    params:
        db=KRAKEN,
    conda:
        "envs/kraken.yml"
    log:
        os.path.join(OUT,"report/kraken_{file}.log"),
    benchmark:
        os.path.join(OUT,"benchmark/kraken_{file}.benchmark"),
    output:
        out=os.path.join(OUT,"Pangenome/Kraken/{file}.out"),
        report=os.path.join(OUT,"Pangenome/Kraken/{file}.report")
    shell:
        """
        kraken2 --db {params.db} --output {output.out} --threads 20 --use-names --report {output.report} --use-mpa-style {input.file} &> {log}
        """

rule kraken_reformat:
    input:
        out=os.path.join(OUT,"Pangenome/Kraken/{file}.out"),
    params:
        sample="{file}",
        outpath=os.path.join(OUT,"Pangenome/Kraken"),
    log:
        os.path.join(OUT,"report/kraken_reformat_{file}.log"),
    benchmark:
        os.path.join(OUT,"benchmark/kraken_reformat_{file}.benchmark"),
    output:
        report=os.path.join(OUT,"Pangenome/Kraken/{file}_reformat.txt"),
    shell:
        "python panpiper/workflow/scripts/kraken_reformat.py {params.sample} {params.outpath} &> {log}"

rule kraken_ag:
    input:
        report=expand(os.path.join(OUT,"Pangenome/Kraken/{file}_reformat.txt"), file=filtered),
    output:
        os.path.join(OUT,"Pangenome/Summary/kraken_ag.txt"),
    shell:
        "awk FNR!=1 {input} > {output}"


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

