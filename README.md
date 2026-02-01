[![Project Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

# Panpiper

**Panpiper** is a comprehensive Snakemake-based bioinformatics pipeline for bacterial isolate analysis. It provides an end-to-end workflow from raw sequencing reads through genome assembly, quality control, pangenome construction, and genome-wide association studies (GWAS).

## Features

- **Modular Design**: Run individual workflows independently or as a complete pipeline
- **Reproducible Analysis**: Snakemake-based workflow with isolated Conda environments
- **Scalable**: Supports local execution and HPC cluster environments
- **Comprehensive Output**: Generates publication-ready results including phylogenetic trees, gene matrices, and GWAS outputs

## Workflows

### Assembly
Assembles paired-end FASTQ files into contigs using [Shovill](https://github.com/tseemann/shovill).

### Quality Control
Filters assemblies based on quality metrics using [CheckM2](https://github.com/chklovski/CheckM2) and validates taxonomy with [FastANI](https://github.com/ParBLiSS/FastANI).

### Pangenome Analysis
Constructs the pangenome using [Panaroo](https://github.com/gtonkinhill/panaroo) with functional annotation via:
- [Bakta](https://github.com/oschwengers/bakta) - Gene prediction and annotation
- [AMRFinderPlus](https://github.com/ncbi/amr) - Antimicrobial resistance gene detection
- [EggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) - Functional annotation

Phylogenetic trees are constructed using [FastTree](http://www.microbesonline.org/fasttree/), [RAxML](https://github.com/stamatak/standard-RAxML), and [IQ-TREE](https://github.com/Cibiv/IQ-TREE).

### Genome-Wide Association Study
Performs genotype-phenotype association analysis using [Pyseer](https://github.com/mgalardini/pyseer) with support for:
- Gene presence/absence analysis
- Structural variant analysis
- Unitig (k-mer) analysis

## Installation

### Prerequisites
- Conda or Mamba
- Python 3.8+

### Setup

```bash
# Create and activate environment
conda create -n panpiper -c bioconda -c conda-forge snakemake mamba
conda activate panpiper

# Install Panpiper
pip install panpiper
```

### Database Configuration

Most databases are downloaded automatically. For custom database locations:

```bash
panpiper ... --checkm2_dir /path/to/checkm2 --eggnog_dir /path/to/eggnog --bakta_dir /path/to/bakta
```

## Usage

### Assembly

Assemble paired-end FASTQ files:

```bash
panpiper -w assembly -o output_dir/ -q fastq_dir/
```

**Input Requirements:**
- Directory containing paired-end FASTQ files (`*_1.fastq.gz`, `*_2.fastq.gz`)

**Output:**
| File | Description |
|------|-------------|
| `Assembly/{sample}/contigs.fa` | Assembled contigs per sample |
| `Assembly/complete_assembly_files.txt` | Samples with successful assemblies |
| `Assembly/incomplete_assembly_files.txt` | Samples with failed assemblies |

### Quality Control

Filter assemblies by quality metrics and taxonomic verification:

```bash
panpiper -w quality -o output_dir/ -a fasta_dir/ -s sample_list.txt -r reference.fasta
```

**Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--ani_cutoff` | 95 | Minimum ANI to reference (%) |
| `--n50` | 5000 | Minimum N50 value |
| `--contig_number` | 1000 | Maximum number of contigs |

**Output:**
| File | Description |
|------|-------------|
| `Quality/sample_list.txt` | Samples passing all criteria |
| `Quality/failed_samples_checkm.csv` | Samples failing CheckM thresholds |
| `Quality/failed_samples_ani.csv` | Samples failing ANI threshold |
| `Quality/quality_report.html` | Visual quality summary |
| `Quality/CheckM/checkm_cat.txt` | Combined CheckM statistics |
| `Quality/FastANI/fastani_summary.txt` | ANI results per sample |

### Pangenome Analysis

Construct and analyze the pangenome:

```bash
panpiper -w pangenome -o output_dir/ -a fasta_dir/ -s sample_list.txt -r reference.fasta
```

**Output:**
| Directory/File | Description |
|----------------|-------------|
| `Pangenome/Bakta/` | Per-sample gene annotations |
| `Pangenome/Panaroo/` | Pangenome matrices and alignments |
| `Pangenome/AMR/` | Antimicrobial resistance reports |
| `Pangenome/Unitig/unitig.pyseer` | Unitig output for GWAS |
| `Pangenome/Summary/core_gene_alignment.aln.iqtree` | Final phylogenetic tree (Newick) |
| `Pangenome/Summary/genes_matrix.txt` | Gene presence/absence matrix |
| `Pangenome/Summary/genes_anno.txt` | Gene annotations |

### Genome-Wide Association Study

Perform association analysis between genotype and phenotype:

```bash
panpiper -w pyseer -o output_dir/ \
    -g gene_presence_absence.Rtab \
    -p struct_presence_absence.Rtab \
    -u unitig.pyseer.gz \
    -t tree.nwk \
    -f phenotypes.txt \
    -r references.txt \
    --pheno_column trait_name
```

**Input Requirements:**
- Gene presence/absence matrix (`.Rtab` from Panaroo)
- Structure presence/absence matrix
- Unitig file (gzipped)
- Phylogenetic tree (Newick format)
- Phenotype file (tab-delimited with sample names)
- Reference file for unitig annotation

**Output:**
| File | Description |
|------|-------------|
| `Pyseer/{phenotype}/gene_analysis.txt` | Gene association results |
| `Pyseer/{phenotype}/significant_genes.txt` | Significant gene hits |
| `Pyseer/{phenotype}/unitig.txt` | Unitig association results |
| `Pyseer/{phenotype}/significant_unitig.txt` | Significant unitig hits |
| `Pyseer/{phenotype}/unitig_gene_hits.txt` | Annotated significant unitigs |
| `Pyseer/{phenotype}/struct_analysis.txt` | Structural variant results |

## Cluster Execution

For HPC cluster execution, provide a Snakemake profile:

```bash
panpiper ... --cluster --profile profile_dir/
```

The profile directory must contain a `config.yaml` file. See [Snakemake profiles documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) for configuration details.

## Visualization

Results can be visualized using the interactive [Panpiper Streamlit App](https://panpiper.streamlit.app/).

## Citation

If you use Panpiper in your research, please cite the individual tools used in the pipeline:
- Shovill, CheckM2, FastANI, Panaroo, Bakta, AMRFinderPlus, EggNOG-mapper, FastTree, RAxML, IQ-TREE, Pyseer

## Acknowledgments

Pipeline architecture inspired by [MAGinator](https://github.com/Russel88/MAGinator).

## License

This project is licensed under the BSD 3-Clause License - see the [LICENSE](LICENSE) file for details.
