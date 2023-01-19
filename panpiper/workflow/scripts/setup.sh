#!/bin/bash

SAMPLES="${*%${!#}}"
OUT_DIR="${!#}" 

# Fsm-lite
#mkdir ../MGWAS/Fsm-lite
#exec 4<"$1"
#echo Start
#while IFS=$' \t\r\n' read -r p; do
#	echo -e "$p ../Assembly/Shovill/$p/$p.fna"
#done < $1 > ../MGWAS/Fsm-lite/fsm-lite.list

# Poppunk
mkdir -p $OUT_DIR/Phylogroups
exec 4<"$1"
echo Building poppunk list
for p in $SAMPLES; do 
	echo $p
	NAME=basename $p .fna
	echo -e "$NAME	$p" >> $OUT_DIR/Phylogroups/poppunk.list
done

# Snippy
#mkdir ../MGWAS/Snippy
#exec 4<"$1"
#echo Start
#while IFS=$' \t\r\n' read -r p; do
#	echo -e "$p	../Assembly/Shovill/$p/$p.fna"
#done < $1 > ../MGWAS/Snippy/snippy.list

mkdir -p $OUT_DIR/Unitig
# Unitig
exec 4<"$1"
echo Building unitig list
for p in $SAMPLES; do 
	echo $p
	echo -e $p >> $OUT_DIR/Unitig/unitig.list
done

#for f in *_filtered_sorted; do 
#	n=$(echo $f | cut -d'_' -f 1)
#	mv "$f" "$n"
#done
