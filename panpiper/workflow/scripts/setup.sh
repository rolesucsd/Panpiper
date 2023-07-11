SAMPLES="${*%${!#}}"
OUT_DIR="${!#}" 

mkdir -p $OUT_DIR/Unitig
# Unitig
exec 4<"$1"
echo Building unitig list
for p in $SAMPLES; do 
	echo $p
	echo -e $p >> $OUT_DIR/Unitig/unitig.list
done