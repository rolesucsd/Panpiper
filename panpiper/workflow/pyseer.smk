# CHECKPOINT #
"""
Authors: Renee Oles
Date Updated: 1/4/2022
Purpose: Genome-wide association study

TODO:
- make sure you have constructed phenotypes as either binary (needs to be 1 or 0) or continuous
- see if there is a way to include multiple phenotypes together
- include a file for covariates (confounding variables) such as the group that the sample is from
- need to make the reference.txt file automatically using prokka files - and sample list output from isoqual
- add power calculations to the results so that the user can see how much power they have to make conclusions 
- need to change mash output so that it has just the sample name and not the full path
"""

configfile: "config/config.yaml"

rule all:
    input:
#        ("full/ref_unitig.plot"),
#        ("full/unitig_pattern_count.txt"),
#        ("full/SNPs.txt"),
        ("growth/COG_analysis.txt"),
#        ("full/SNPs_lineage1.txt"),
#        ("full/unitig_qqplot.png"),
        ("growth/unitig_gene_hits.txt"),
        ("growth/struct_analysis.txt"),
        ("growth/unitig_gene_hits_enet.txt")

# We will use the distances from the core genome phylogeny, which has been midpointed rooted:
# TO DO: Has the tree been midpoint rooted though????? 
rule pyseer_SNP_to_matrix:
    input:
        "full/core.tre"
    conda:
        "envs/pyseer.yml"
    output:
        "full/phylogeny_similarity.tsv"
    shell:
        "python scripts/phylogeny_distance.py --lmm {input} > {output}"

rule pyseer_SNP_to_dist:
    input:
        "full/core.tre"
    conda:
        "envs/pyseer.yml"
    output:
        "full/phylogeny_distance.tsv"
    shell:
        "python scripts/phylogeny_distance.py {input} > {output}"


# TODO: May want to include covariates 
rule pyseer_COG_analysis:
    input:
        roary="full/gene_presence_absence.Rtab",
        phylo="full/phylogeny_similarity.tsv"
    params:
        pheno_col=config["pheno"],
        pheno=config["meta"],
        dimensions=2
#        cov_col=config["cov"]
    conda:
        "envs/pyseer.yml"  
    output:
        COG="growth/COG_analysis.txt"
    shell:
        "pyseer --lmm --phenotypes {params.pheno} --phenotype-column {params.pheno_col} --pres {input.roary} --similarity {input.phylo} --cpu 8 > {output.COG}"
       
rule pyseer_structure_analysis:
    input:
        roary="full/struct_presence_absence.Rtab",
#        mash="full/fasta.tsv",
        phylo="full/phylogeny_similarity.tsv"
    params:
        pheno_col=config["pheno"],
        pheno=config["meta"],
        dimensions=2
#        cov_col=config["cov"]
    conda:
        "envs/pyseer.yml"  
    output:
#        patterns="cog_patterns.txt",
        COG="growth/struct_analysis.txt"
    shell:
        "pyseer --lmm --phenotypes {params.pheno} --phenotype-column {params.pheno_col} --pres {input.roary} --similarity {input.phylo}  > {output.COG}"
        #--distances {input.mash} --covariates {input.pheno} --use-covariates {params.cov_col} --save-m ../Pyseer/mash_mds --max-dimensions {params.dimensions}
#        "pyseer --phenotypes {input.pheno} --phenotype-column {params.pheno_col} --pres {input.roary} --no-distances --output-patterns {output.patterns} --cpu 8 > {output.COG}"

rule pyseer_unitig:
    input:
        pheno=config["meta"],
        phylo="full/phylogeny_similarity.tsv",
        unitig="full/unitig.pyseer.gz",
        phylo_dist="full/phylogeny_distance.tsv"
    params:
        pheno_col=config["pheno"],
        cov_col=config["cov"]
    conda:
        "envs/pyseer.yml"
    output:
        patterns="growth/unitig_patterns.txt",
        unitig="growth/unitig.txt"
    shell:
#        pyseer s --print-samples --phenotypes {input.pheno} --phenotype-column {params.pheno_col} --covariates {input.pheno} --use-covariates {params.cov_col} --kmers {input.fsm} --no-distance --output-patterns {output.patterns} --cpu 8 > {output.kmers}
#        gzip {input.unitig}
        """
        pyseer --lmm --print-samples  --similarity {input.phylo} --phenotypes {input.pheno} --phenotype-column {params.pheno_col} --kmers {input.unitig} --output-patterns {output.patterns} --cpu 8 --distances {input.phylo_dist} --lineage > {output.unitig}
        """

rule unitig_qqplot:
    input:
        "growth/unitig.txt"
    conda:
        "envs/pyseer.yml"
    output:
        "growth/unitig_qqplot.png"
    shell:
        "python scripts/qq_plot.py {input} --output {output}"

rule pyseer_patterns_unitig:
    input:
        "growth/unitig_patterns.txt"
    conda:
        "envs/pyseer.yml"
    output:
        "growth/unitig_pattern_count.txt"
    shell:
        """
        python scripts/count_patterns.py {input} > {output}
        """

# TO DO: Need to make a python script to do the same thing as the command below 
rule pyseer_filter_unitig:
    input:
        n= "growth/unitig_pattern_count.txt",
        nxt="growth/unitig.txt"
    conda:
        "envs/r.yml"
    output:
        "growth/significant_unitig.txt"
    shell:
        """
        Rscript scripts/filter_kmers.R {input.nxt} {output}
        """

# Using bwa mem to map significant kmers to reference provided
rule pyseer_map_to_ref_unitig:
    input:
        unitig="growth/significant_unitig.txt",
        inte="pan_reference/pan_genome_reference.fna.bwt"
    params:
        ref="growth/pan_genome_reference.fa"
    conda:
        "envs/pyseer.yml"
    output:
        "growth/ref_unitig.plot"
    shell:
        "phandango_mapper {input.unitig} {params.ref} {output}"

rule annotate_unitig:
    input:
        unitig="growth/significant_unitig.txt",
        bwt1="metadata/d1.fna.bwt",
        bwt2="metadata/d2.fna.bwt"
    params:
        "metadata/reference.txt"
    conda:
        "envs/pyseer.yml"
    output:
        "growth/unitig_annotation.txt"
    shell:
        "annotate_hits_pyseer {input.unitig} {params} {output}"

rule unitig_gene:
    input:
        "growth/unitig_annotation.txt"
    conda:
        "envs/pyseer.yml"
    output:
        "growth/unitig_gene_hits.txt"
    shell:
        "python scripts/summarise_annotations.py {input} > {output}"

rule pyseer_elastic_net:
    input:
        pheno=config["meta"],
        unitig="full/unitig.pyseer.gz"
    params:
        pheno_col=config["pheno"],
        cov_col=config["cov"]
    conda:
        "envs/pyseer.yml"
    output:
        unitig="growth/unitig_enet.txt"
    shell:
        "pyseer --wg enet --cpu 8 --phenotypes {input.pheno} --phenotype-column {params.pheno_col} --kmers {input.unitig}"

rule pyseer_filter_unitig_elastic_net:
    input:
        nxt="growth/unitig_enet.txt"
    conda:
        "envs/r.yml"
    output:
        "growth/significant_unitig_enet.txt"
    shell:
        """
        Rscript scripts/filter_kmers.R {input.nxt} {output}
        """

rule annotate_unitig_elastic_net:
    input:
        unitig="growth/significant_unitig_enet.txt",
        bwt1="metadata/d1.fna.bwt",
        bwt2="metadata/d2.fna.bwt"
    params:
        "metadata/reference.txt"
    conda:
        "envs/pyseer.yml"
    output:
        "growth/unitig_annotation_enet.txt"
    shell:
        "annotate_hits_pyseer {input.unitig} {params} {output}"

rule unitig_gene_elastic_net:
    input:
        "growth/unitig_annotation_enet.txt"
    conda:
        "envs/pyseer.yml"
    output:
        "growth/unitig_gene_hits_enet.txt"
    shell:
        "python scripts/summarise_annotations.py {input} > {output}"


