[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

# Package

Keyword description: 

* Accomplishes 1
    * Explanation
* Accomplishes 2
    * Explanation

## Installation

Create package and install dependencies 

```sh
conda create -n package -c bioconda -c conda-forge snakemake mamba
conda activate package
pip install package
```

Download any necessary databases (perhaps kraken)
```sh
code
```

## Usage

Input: 

* subdirectory of the package directory intitled "resources"
    * directory fastq for assembly workflow
    * directory fasta for quality workflow
    * directory ref for quality workflow
    * directory fasta for pangenome workflow 

Run Package:

* Package can be run at three starting points: 

```sh
package -q fastq directory -a fasta directory -s sample_list -r reference file with path -o output directory "
```


### Run on a compute cluster
MAGinator can run on compute clusters using qsub (torque), sbatch (slurm), or drmaa structures. The --cluster argument toggles the type of compute cluster infrastructure. The --cluster_info argument toggles the information given to the submission command, and it has to contain the following keywords {cores}, {memory}, {runtime}, which are used to forward resource information to the cluster.

A qsub MAGinator can for example be run with the following command (... indicates required arguments, see above):
```sh
maginator ... --cluster qsub --cluster_info "-l nodes=1:ppn={cores}:thinnode,mem={memory}gb,walltime={runtime}"
```

## Workflow


## Output

* folder/
    * file - description
    
