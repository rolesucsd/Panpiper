#!/bin/bash

if [ -f $1/complete_assembly_files.txt ]; then
	rm $1/complete_assembly_files.txt
	touch $1/complete_assembly_files.txt
else
	touch $1/complete_assembly_files.txt
fi

if [ -f $1/incomplete_assembly_files.txt ]; then
	rm $1/incomplete_assembly_files.txt
	touch $1/incomplete_assembly_files.txt
else
	touch $1/incomplete_assembly_files.txt
fi

CONTIGS="contigs.fa"
for d in $1/*/
do
	if [ -s $d$CONTIGS ]; then
		echo "$(basename -- $d)" >> $1/complete_assembly_files.txt
	else
		echo "$(basename -- $d)" >> $1/incomplete_assembly_files.txt
	fi
done
