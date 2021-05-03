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
trimmomatic PE -threads $THREADS -phred33 -summary $STATS_TRIM/${i}.summary \
     ${i}$EXTENSION_PE1 ${i}$EXTENSION_PE2 \
     $TRIM/${i}.R1_paired.fastq.gz $TRIM/${i}.R1_unpaired.fastq.gz \
     $TRIM/${i}.R2_paired.fastq.gz $TRIM/${i}.R2_unpaired.fastq.gz \
     ILLUMINACLIP:"$ADAPTER_TRIMMOMATIC" \
     LEADING:15 TRAILING:15 SLIDINGWINDOW:5:20 AVGQUAL:20 MINLEN:$MIN_READ_LENGTH CROP:$MAX_READ_LENGTH


END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ TRIMMOMATIC SAMPLE $i TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "#@##################################################################"

echo
echo
echo
