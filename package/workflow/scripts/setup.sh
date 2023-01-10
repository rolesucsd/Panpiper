#!/bin/bash

# Fsm-lite
mkdir ../MGWAS/Fsm-lite
exec 4<"$1"
echo Start
while IFS=$' \t\r\n' read -r p; do
	echo -e "$p ../Assembly/Shovill/$p/$p.fna"
done < $1 > ../MGWAS/Fsm-lite/fsm-lite.list

# Poppunk
mkdir ../MGWAS/Phylogroups
exec 4<"$1"
echo Start
while IFS=$' \t\r\n' read -r p; do
	echo -e "$p	../Assembly/Shovill/$p/$p.fna"
done < $1 > ../MGWAS/Phylogroups/poppunk.list

# Snippy
#mkdir ../MGWAS/Snippy
#exec 4<"$1"
#echo Start
#while IFS=$' \t\r\n' read -r p; do
#	echo -e "$p	../Assembly/Shovill/$p/$p.fna"
#done < $1 > ../MGWAS/Snippy/snippy.list

mkdir ../MGWAS/Unitig
# Unitig
exec 4<"$1"
echo Start
while IFS=$' \t\r\n' read -r p; do
	echo -e "../Assembly/Shovill/$p/$p.fna"
done < $1 > ../MGWAS/Unitig/unitig.list

#for f in *_filtered_sorted; do 
#	n=$(echo $f | cut -d'_' -f 1)
#	mv "$f" "$n"
#done
