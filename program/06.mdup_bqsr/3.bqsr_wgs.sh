#!/bin/bash


echo "#@##################################################################"
echo "#@ BASE QUALITY SCORE RECALIBRATION: SAMPLE ${i}"
echo
D1=`date "+%D    %T"`
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

gatk --java-options '-Djava.io.tmpdir=$TMP -Xmx20g' BaseRecalibrator \
   -R $REFERENCE \
   -I $ALIGN_MDUP/${i}_trim_"$ALIGNER"_sorted_RG_mdup.bam \
   -known-sites $VARIANTS \
   -O $ALIGN_BQSR/${i}_recalibration.table \
   --tmp-dir $TMP

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ BQSR SAMPLE ${i} TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "#@##################################################################"

echo
echo
echo

echo "#@##################################################################"
echo "#@ APPLY BQSR: SAMPLE ${i}"
echo
D2=`date "+%D    %T"`
echo "#@ Date and Time: $D2"
echo
       #################################


START2=$(date +%s)

gatk --java-options '-Djava.io.tmpdir=$TMP' ApplyBQSR \
   -R $REFERENCE \
   -I $ALIGN_MDUP/${i}_trim_"$ALIGNER"_sorted_RG_mdup.bam \
   --bqsr-recal-file $ALIGN_BQSR/${i}_recalibration.table \
   -O $ALIGN_BQSR/${i}_trim_"$ALIGNER"_sorted_RG_mdup_bqsr.bam \
   --tmp-dir $TMP
wait

rm -r $ALIGN_BQSR/${i}_recalibration.table 2> /dev/null
rm -r $ALIGN_BQSR/${i}_trim_"$ALIGNER"_sorted_RG_mdup_bqsr.bam.sbi 2> /dev/null
wait

END2=$(date +%s)
DIFF2=$(( $END2 - $START2 ))

echo
echo "#@ APPLY BQSR SAMPLE ${i} TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF2/3600)) $(($DIFF2%3600/60)) $(($DIFF2%60)))"
echo
echo "#@##################################################################"

echo
