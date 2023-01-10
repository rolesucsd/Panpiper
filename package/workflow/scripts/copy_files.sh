#!/bin/bash

mkdir resources/fasta
while IFS= read -r line
do
    cp "results/Assembly/{$line}/{$line.fna}" "resources/fasta/{$line}.fna"
done < results/Quality/sample_list.txt
