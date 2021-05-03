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

if [[ $TYPE == "se" ]]; then
   featureCounts -T $THREADS -t $FEATURE -a $ANNOTATION -o $GENECOUNT_DATA/${i}_counts_"$FEATURE".txt $ALIGN/${i}_trim_"$ALIGNER"_sorted_RG.bam
else
   featureCounts -p -T $THREADS -t $FEATURE -a $ANNOTATION -o $GENECOUNT_DATA/${i}_counts_"$FEATURE".txt $ALIGN/${i}_trim_"$ALIGNER"_sorted_RG.bam
fi
wait

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ $ALIGNER $GENE_COUNT $FEATURE SAMPLE $i TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "#@##################################################################"

