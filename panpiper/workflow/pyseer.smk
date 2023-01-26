# ----------------------------------------------------------------------------
# Copyright (c) 2022--, Panpiper development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""
Authors: Renee Oles
Purpose: Genome-wide association study
TODO:
- make sure you have constructed phenotypes as either binary (needs to be 1 or 0) or continuous
- see if there is a way to include multiple phenotypes together
- include a file for covariates (confounding variables) such as the group that the sample is from
- need to make the reference.txt file automatically using prokka files - and sample list output from isoqual
- add power calculations to the results so that the user can see how much power they have to make conclusions 
- need to change mash output so that it has just the sample name and not the full path
"""
import os
import subprocess
import math
import distutils.util

OUT = config['out']
COG = config['cog']
STRUCTURE = config['structure']
UNITIG = config['unitig']
SNP = config['snp']
PARAMS = config['params']
REF = config['ref']

with open(PARAMS, 'r') as fh:   
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

rule all:
    input:
        os.path.join(OUT,"Pyseer/COG_analysis.txt"),
        os.path.join(OUT,"Pyseer/unitig_gene_hits.txt"),
        os.path.join(OUT,"Pyseer/struct_analysis.txt"),
        os.path.join(OUT,"Pyseer/unitig_gene_hits_enet.txt"),

# We will use the distances from the core genome phylogeny, which has been midpointed rooted:
# TODO: check root of tree
rule pyseer_SNP_to_matrix:
    input:
        SNP
    conda:
        "envs/pyseer.yml"
    output:
        os.path.join(OUT,"Pyseer/phylogeny_similarity.tsv"),
    shell:
        "python panpiper/workflow/scripts/phylogeny_distance.py --lmm {input} > {output}"

# TODO: May want to include covariates 
rule pyseer_COG_analysis:
    input:
        cog=COG
        phylo=os.path.join(OUT,"Pyseer/phylogeny_similarity.tsv"),
    params:
        pheno_col=param_dict["pheno"],
        pheno=param_dict["meta"],
    conda:
        "envs/pyseer.yml"  
    output:
        os.path.join(OUT,"Pyseer/COG_analysis.txt"),
    shell:
        "pyseer --lmm --phenotypes {params.pheno} --phenotype-column {params.pheno_col} --pres {input.roary} --similarity {input.phylo} --cpu 8 > {output}"
       
rule pyseer_structure_analysis:
    input:
        roary=STRUCTURE,
        phylo=os.path.join(OUT,"Pyseer/phylogeny_similarity.tsv"),
    params:
        pheno_col=param_dict["pheno_column"],
        pheno=param_dict["pheno_file"],
    conda:
        "envs/pyseer.yml"  
    output:
        COG=os.path.join(OUT,"Pyseer/struct_analysis.txt")
    shell:
        "pyseer --lmm --phenotypes {params.pheno} --phenotype-column {params.pheno_col} --pres {input.roary} --similarity {input.phylo}  > {output.COG}"

rule pyseer_unitig:
    input:
        phylo=os.path.join(OUT,"Pyseer/phylogeny_similarity.tsv"),
        unitig=UNITIG,
    params:
        pheno_col=param_dict["pheno_column"],
        pheno=param_dict["pheno_file"],
    conda:
        "envs/pyseer.yml"
    output:
        patterns=os.path.join(OUT,"Pyseer/unitig_patterns.txt"),
        unitig=os.path.join(OUT,"Pyseer/unitig.txt"),
    shell:
        "pyseer --lmm --print-samples  --similarity {input.phylo} --phenotypes {params.pheno} --phenotype-column {params.pheno_col} --kmers {input.unitig} --output-patterns {output.patterns} --cpu 8 --lineage > {output.unitig}"

rule pyseer_patterns_unitig:
    input:
        os.path.join(OUT,"Pyseer/unitig_patterns.txt"),
    conda:
        "envs/pyseer.yml"
    output:
        os.path.join(OUT,"Pyseer/unitig_pattern_count.txt"),
    shell:
        "python panpiper/workflow/scripts/count_patterns.py {input} > {output}"


rule pyseer_filter_unitig:
    input:
        n=os.path.join(OUT,"Pyseer/unitig_pattern_count.txt"),
        nxt=os.path.join(OUT,"Pyseer/unitig.txt"),
    conda:
        "envs/r.yml"
    output:
        os.path.join(OUT,"Pyseer/significant_unitig.txt"),
    shell:
        "Rscript panpiper/workflow/scripts/filter_kmers.R {input.nxt} {output}"


rule annotate_unitig:
    input:
        unitig=os.path.join(OUT,"Pyseer/significant_unitig.txt")
    params:
        REF
    conda:
        "envs/pyseer.yml"
    output:
        unitig=os.path.join(OUT,"Pyseer/unitig_annotation.txt"),
    shell:
        "annotate_hits_pyseer {input.unitig} {params} {output}"

rule unitig_gene:
    input:
        os.path.join(OUT,"Pyseer/unitig_annotation.txt"),
    conda:
        "envs/pyseer.yml"
    output:
       os.path.join(OUT,"Pyseer/unitig_gene_hits.txt"),
    shell:
        "python panpiper/workflow/scripts/summarise_annotations.py {input} > {output}"

rule pyseer_elastic_net:
    input:
        unitig=os.path.join(OUT,"Pyseer/unitig.pyseer.gz"),
    params:
        pheno_col=param_dict["pheno_column"],
        pheno=param_dict["pheno_file"],
    conda:
        "envs/pyseer.yml"
    output:
        unitig=os.path.join(OUT,"Pyseer/unitig_enet.txt"),
    shell:
        "pyseer --wg enet --cpu 8 --phenotypes {params.pheno} --phenotype-column {params.pheno_col} --kmers {input.unitig}"

rule pyseer_filter_unitig_elastic_net:
    input:
        nxt=os.path.join(OUT,"Pyseer/unitig_enet.txt"),
    conda:
        "envs/r.yml"
    output:
        unitig=os.path.join(OUT,"Pyseer/significant_unitig_enet.txt"),
    shell:
        """
        Rscript panpiper/workflow/scripts/filter_kmers.R {input.nxt} {output}
        """

rule annotate_unitig_elastic_net:
    input:
        unitig=os.path.join(OUT,"Pyseer/significant_unitig_enet.txt"),
    params:
        REF
    conda:
        "envs/pyseer.yml"
    output:
        os.path.join(OUT,"Pyseer/unitig_annotation_enet.txt"),
    shell:
        "annotate_hits_pyseer {input.unitig} {params} {output}"

rule unitig_gene_elastic_net:
    input:
       os.path.join(OUT,"Pyseer/unitig_annotation_enet.txt"),
    conda:
        "envs/pyseer.yml"
    output:
        os.path.join(OUT,"Pyseer/unitig_gene_hits_enet.txt"),
    shell:
        "python panpiper/workflow/scripts/summarise_annotations.py {input} > {output}"

