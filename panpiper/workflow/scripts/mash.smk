# Mash- computes approximate distance between two samples QUICKLY 
# locality-sensitive hashing (called MinHash sketches which have been used in lots of areas ie website comparisons)
# Time: 30 minutes for 571 samples
rule mash_fasta:
    input:
        file=expand(os.path.join(FASTA, "{file}/{file}.fna"), file=filtered),
    params:
        out=os.path.join(OUT, "Pangenome/Summary/mash"),
    conda:
        "envs/mash.yml",
    output:
        file=os.path.join(OUT, "Pangenome/Summary/mash.tsv"),
    shell:
        """
        mash sketch -s 10000 -o {params.out} {input}
        mash dist {params.out}.msh {params.out}.msh -t > {output}
        """
#        sed -i -e 's/.fasta//g' {output}
#       sed -i -e "s/${params.pref}//g" {output}"
