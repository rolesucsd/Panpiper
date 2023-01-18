#!/bin/bash

SAMPLES="${*%${!#}}"
OUT_FILE="${!#}" 

if [ -f $OUT_FILE ]; then
    rm $OUT_FILE
    touch $OUT_FILE
else
    touch $OUT_FILE
fi

for i in $SAMPLES; do 
    echo "hi"
    echo "$i"  >> $OUT_FILE
 done