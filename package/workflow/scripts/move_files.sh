#!/bin/bash

for i in $1
do
	if [ -s $1$i/contigs.fa ]; then
		mv $1$i/contigs.fa $1complete/$i
	fi
done

ls $1complete/* > $1complete_assembly_files.txt