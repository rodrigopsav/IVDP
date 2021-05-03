#!/bin/bash


echo "#@##################################################################"
echo "#@ FASTQC RAW FILES: SAMPLE ${i}"
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

cd $INPUT_DIR
fastqc -o $STATS_QC1 -t $THREADS -dir $TMP ${i}$EXTENSION_PE1 ${i}$EXTENSION_PE2

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ FASTQC SAMPLE $i TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "#@##################################################################"

echo
echo
echo