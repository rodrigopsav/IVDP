#!/bin/bash


echo "#@##################################################################"
echo "#@ SAMPLE: $i --- ADD RG GROUP START"
echo "#@ ADD OR REPLACE READ GROUPS"
echo
D1=`date "+%D    %T"`
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

cd $ALIGN
picard -Xmx20g AddOrReplaceReadGroups \
   I=${i}_trim_"$ALIGNER"_sorted.bam \
   O=${i}_trim_"$ALIGNER"_sorted_RG.bam \
   RGID=$ID_PU \
   RGLB=lib1 \
   RGPL=illumina \
   RGPU=$ID_PU \
   RGSM=$SUBJECT_ID \
   VALIDATION_STRINGENCY=LENIENT

cd $ALIGN
samtools index -@ $THREADS -b ${i}_trim_"$ALIGNER"_sorted_RG.bam
wait

rm $ALIGN/${i}_trim_"$ALIGNER"_sorted.bam
rm $ALIGN/${i}_trim_"$ALIGNER"_sorted.bam.bai

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ ADD RG SAMPLE $i TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "#@##################################################################"

echo
