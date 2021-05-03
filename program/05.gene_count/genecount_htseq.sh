#!/bin/bash


echo "#@##################################################################"
echo "#@ GENE COUNT"
echo "#@ $ALIGNER $GENE_COUNT $FEATURE: SAMPLE $i"
echo
D1=`date "+%D    %T"`
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

htseq-count -f bam -r pos -t $FEATURE -m union -s no $ALIGN/${i}_trim_"$ALIGNER"_sorted_RG.bam $ANNOTATION > $GENECOUNT_DATA/${i}_counts_"$FEATURE".txt
wait

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ $ALIGNER $GENE_COUNT $FEATURE SAMPLE $i TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "#@##################################################################"

echo
