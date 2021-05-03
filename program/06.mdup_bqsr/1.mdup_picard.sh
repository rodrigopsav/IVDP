#!/bin/bash

############ PATH FOLDER ###########
mkdir -p $STATS_ALIGN_CORRECTION 2>&1 > /dev/null
####################################

echo "#@##################################################################"
echo "#@ MARK DUPLICATED READS: SAMPLE $i"
echo
D1=`date "+%D    %T"`
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

picard -Xmx20g MarkDuplicates \
   I=$ALIGN/"$i"_trim_"$ALIGNER"_sorted_RG.bam \
   O=$ALIGN_MDUP/"$i"_trim_"$ALIGNER"_sorted_RG_mdup.bam \
   M=$STATS_ALIGN/"$i"_mdup_picard_metrics.txt
wait

samtools index -@ $THREADS -b $ALIGN_MDUP/"$i"_trim_"$ALIGNER"_sorted_RG_mdup.bam

#rm -v $ALIGN/"$i"_trim_"$ALIGNER"_sorted_RG.bam
#rm -v $ALIGN/"$i"_trim_"$ALIGNER"_sorted_RG.bam.bai


END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ MDUP SAMPLE $i TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "#@##################################################################"

echo

