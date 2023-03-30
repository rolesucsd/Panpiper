[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

# Panpiper

This package conducts bacterial isolate analysis for one species.  

* Assembly
    * Using [shovill](https://github.com/tseemann/shovill) the paired fastq files are assembled and sorted into ones which pass and fail assembly
* Quality control
    * The assemblies are analyzed with [CheckM](https://github.com/Ecogenomics/CheckM) and [FastANI](https://github.com/ParBLiSS/FastANI)
    * Based on user-defined thresholds the samples are sorted into ones which pass or don't pass requirements
* Pangenome analysis
    * A pangenome is created with the assemblies using [Panaroo](https://github.com/gtonkinhill/panaroo) and annotated with [AMRFinderPlus](https://github.com/ncbi/amr), [Bakta](https://github.com/oschwengers/bakta), [EggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper), and [Kraken2](https://github.com/DerrickWood/kraken2)
    * The core genome alignment is used to create a phylogenetic tree with [FastTree](http://www.microbesonline.org/fasttree/), [RAxML](https://github.com/stamatak/standard-RAxML), and [IQ-TREE](https://github.com/Cibiv/IQ-TREE)
    * The samples are divided into phylogroups with [PopPUNK](https://github.com/bacpop/PopPUNK)
* Genome-wide association study
    * The pangenome is used with a continuous or binary phenotype to conduct a genome-wide association study with [Pyseer](https://github.com/mgalardini/pyseer)
    * The wrapper scripts used to analyze and filter the data are predominantly taken from the pyseer package.

Note: Be sure to cite all packages used in the pipeline. 

Credit: This package structure was inspired by the snakemake workflow [MAGinator](https://github.com/Russel88/MAGinator)

![Workflow Diagram](workflow.png)

## Installation

Create package and install dependencies 

```sh
conda create -n panpiper -c bioconda -c conda-forge snakemake mamba
conda activate panpiper
pip install panpiper
```

Download databases

* The databses necessary will be downloaded automatically except in the case of Kraken2 and EggNOG-mapper. These databases have already been downloaded and are maintained in the package structure. In the case that you would like to use your own database, these databases can be downloaded into user-defined directories and can be referenced in the arguments as follows

```sh
panpiper ... --kraken_dir {directory} --eggnog_dir {directory} --bakta_dir {directory}
```

## Workflow

Note: The first time these functions are run, they will take a couple minutes before the workflow starts running because the conda dependendies need to be downloaded before start.

### Assembly: 
The fastq files should all be in a single directory. They need to be paired-end. 

```sh
panpiper -w assembly -o {output directory} -q {fastq directory} 
```

Output:  
* Assembly/  
    * incomplete_assembly_files.txt - samples for which assembly failed
    *  complete_assembly_files.txt - samples for which assemblies passed 
    * {sample}/  
        * contigs.fa - edited assembly 

### Quality control: 
The fasta file directory should contain a subdirectory for each sample where the subdirectory name is the sample of the sample. The fasta file within these directories should be formatted as contigs.fa. This is the standard output from Shovill. The sample list should have one sample name per line which corresponds to the directory names. The reference fasta file should be the representative strain of the species which the samples will be compared to using average nucleotide identity.

Additional parameters that may be included:
* --ani_cutoff : percent identity to reference cutoff (default 95)
* --n50 : N50 cutoff (default 5000)
* --contig_number : number of contigs cutoff (default 1000)

```sh
panpiper -w quality -o {ouput directory} -a {fasta directory} -s  {sample list} -r {reference fasta file}
```

Output:  
* Quality/  
    * sample_list.txt - samples which passed all requirements
    * failed_samples_checkm.csv  - samples which did not meet checkM thresholds
    * failed_samples_ani.csv  - samples which did not meet ANI threshold
    * quality_report.html - html to view the summary of quality of assemblies
    * CheckM/  
        * checkm_stats.txt - statistics summary of assemblies
        * checkm_log.txt  - summary of contamination, completeness, and strain heterogeneity
        * checkm_log_filter.txt - reformat and filter of checkm_log.txt
    * FastANI/
        * fastani_summary.txt - lists percent identity of each sample to the reference
        * fastani_reformat.csv - reformat of the same file 
    * Assembly_filter/  
        * {sample}/ - includes the filtered fasta files removing contigs less than 500bp and all CheckM information

### Pangenome: 
The fasta file directory should contain a subdirectory for each sample where the subdirectory name is the sample of the sample. The fasta file within these directories should be formatted as contigs.fa. This is the standard output from Shovill. The sample list should have one sample name per line which corresponds to the directory names. The reference fasta file should be the representative strain of the species which the samples will be compared to using average nucleotide identity.

```sh
panpiper -w pangenome  -o {ouput directory} -a {fasta directory} -s  {sample list} -r {reference fasta file}
```

Output:  
* Pangenome/  
    * Bakta/  
        * {sample}/ - contains Bakta output for each sample
    * Unitig/  
        * unitig.pyseer.gz - the full unitig summary of all samples
    * Phylogeny/ - contains the intermediate trees 
    * Panaroo/ - contains the full pangenome summary as well as Bakta annotation for the pangenome
    * AMR/ - contains the full AMR report for each sample
    * Kraken/ - contains the full Kraken report for each sample
    * Phylogroups/ - contains the full division process for phylogroups
    * Summary/  
        * AMR.txt - the summary of AMR hits
        * amr.png - visualization of AMR as a heatmap
        * amr_wide.txt - reformat of AMR hits from long to wide
        * core_gene_alignment.aln.iqtree - the final nwk tree
        * db_clusters.csv - the final phylogroup clusters
        * krakn_ag.txt - the full kraken summary
        * mash.tsv - the full mash distance matrix 
        * pan_genome_reference.faa - all the proteins in the pangenome
        * pan_genome_reference.tsv - the list of all proteins in the pangenome
        * Summary.emapper.annotations - Eggnog-Mapper annotation

### Genome-wide association study: 
The gene and structure presence/absence files should be the result of Panaroo or Roary - from Panaroo this is the .RTab file. The unitig file is a result of unitig-caller; this file should be gzipped. The tree file can be any newick tree - if this pipeline has been followed to this point, we would recommend the tree file from iqtree. Finally, the reference file is a tab-delimited list of files to be used for unitig annotation. The format is file.fna file.gff2 {draft, ref}. 

```sh
panpiper -w pyseer  -o {ouput directory} -g {gene presence absence file} -p {structure presence absence file} -u {unitig file} -t {tree file from iqtree} -r {reference file}
```

Output:   
* Pyseer/  
    * {variable_name}/  
        * unitig_pattern_count.txt - a summary of tests conducted for p-value correction
        * unitig_gene_hits.txt - the genes which significant unitig hits match to
        * unitig.txt - the association study results for each unitig
        * significant_unitig.txt - the unitig.txt file filtered for significance 
        * struct_analysis.txt - the associtaiton study results for each structural variant
        * significant_structure.txt - the structural variant results filtered for significance
        * gene_analysis.txt - the associtaiton study results for each gene 
        * significant_genes.txt - the gene results filtered for significance
        * phylogeny_similarity.tsv - a matrix of similarity to use as a covariate to correct for population structure in the model

### Run on a compute cluster
To run on a compute cluster, a profile folder must be specified which should contain a config.yaml file. Details on what to include in the config file can be found [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles).

```sh
panpiper ... --cluster --profile {profile folder}
```


### Potential errors
1. Kraken - "Error loading hash table" - you may need to request more memory 
2. Kraken - database download errors - I suggest downloading the database manually as it can be very glitchy 

