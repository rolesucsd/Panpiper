
# Setup databases to identify viruses
rule virus_setup:
    input:
        expand(os.path.join(FASTA, "{file}/{file}.fna"), file=filtered),
    conda:
        "envs/virsorter.yml",
    log:
        os.path.join(OUT,"report/virsorter_setup.log"),
    output:
        directory("databases/virsorter"),
    shell:
        """
        mkdir {output}
        virsorter setup -d {output} -j 4 &> {log}
        checkv download_database {output} &> {log}
        DRAM-setup.py prepare_databases --skip_uniref --output_dir {output} &> {log}
        """


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
        "virsorter run --keep-original-seq -i {input.fasta} -w {params.outdir} --include-groups dsDNAphage,ssDNA --min-length 5000 --min-score 0.5 -j 28 all &> {log}"


rule checkV:
    input:
        virsorter=os.path.join(OUT,"Pangenome/Virsorter/VS2_1/{file}/final-viral-combined.fa"),
    params:
        outdir=os.path.join(OUT,"Pangenome/Virsorter/Checkv/{file}"),
    conda:
        "envs/virsorter.yml"
    log:
        os.path.join(OUT,"report/checkv_{file}.log"),
    output:
        os.path.join(OUT,"Pangenome/Virsorter/Checkv/{file}/viruses.fna"),
    shell:
        "checkv end_to_end {input} {params.outdir} -t 28 -d databases/virsorter &> {log}"


rule virsorter_run2:
    input:
        os.path.join(OUT,"Pangenome/Virsorter/Checkv/{file}/viruses.fna"),
    params:
        prov=os.path.join(OUT,"Pangenome/Virsorter/Checkv/{file}/proviruses.fna"),
        combined=os.path.join(OUT,"Pangenome/Virsorter/Checkv/{file}/combined.fna"),
        outdir=os.path.join(OUT,"Pangenome/Virsorter/VS2_2/{file}"),
    conda:
        "envs/virsorter.yml"
    log:
        os.path.join(OUT,"report/virsorter_run2_{file}.log"),
    output:
        i=os.path.join(OUT,"Pangenome/Virsorter/VS2_2/{file}/final-viral-combined-for-dramv.fa"),
        v=os.path.join(OUT,"Pangenome/Virsorter/VS2_2/{file}/viral-affi-contigs-for-dramv.tab"),
    shell:
        """
        cat {input} {params.prov} > {params.combined} 
        virsorter run --seqname-suffix-off --viral-gene-enrich-off --provirus-off --prep-for-dramv -i {params.combined} -w {params.outdir} --include-groups dsDNAphage,ssDNA --min-length 5000 --min-score 0.5 -j 28 all &> {log}
        """


rule dramv:
    input:
        i=os.path.join(OUT,"Pangenome/Virsorter/VS2_2/{file}/final-viral-combined-for-dramv.fa"),
        v=os.path.join(OUT,"Pangenome/Virsorter/VS2_2/{file}/viral-affi-contigs-for-dramv.tab"),
    params:
        outint=os.path.join(OUT,"Pangenome/Virsorter/DRAMv/{file}/annotations.tsv"),
        outdir=os.path.join(OUT,"Pangenome/Virsorter/DRAMv/{file}"),
    conda:
        "envs/virsorter.yml"
    log:
        os.path.join(OUT,"report/dramv_{file}.log"),
    output:
        os.path.join(OUT,"Pangenome/Virsorter/DRAMv/{file}/annotations.tsv"),
    shell:
        """
        DRAM-v.py annotate -i {input} -v {input} -o {params} --skip_trnascan --threads 28 --min_contig_size 1000 &> {log}
        DRAM-v.py distill -i {params.outint} -o {params.outdir} &> {log}
        """

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



