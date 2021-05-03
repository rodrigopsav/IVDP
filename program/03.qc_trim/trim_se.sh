#!/bin/bash


echo "#@##################################################################"
echo "#@ TRIMMOMATIC SAMPLE ${i}"
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

cd $INPUT_DIR
trimmomatic SE -threads $THREADS -phred33 -summary $STATS_TRIM/${i}.summary \
	${i}$EXTENSION_SE \
	$TRIM/${i}_trim.fastq.gz \
	ILLUMINACLIP:"$ADAPTER_TRIMMOMATIC" \
	LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:$MIN_READ_LENGTH CROP:$MAX_READ_LENGTH

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ TRIMMOMATIC SAMPLE $i TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "#@##################################################################"

echo
echo
echo
