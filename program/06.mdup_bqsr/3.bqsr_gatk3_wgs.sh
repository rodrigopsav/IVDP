#!/bin/bash

#eval "$(conda shell.bash hook)"
#conda activate ivdp
#conda activate --stack gatk3

FileName=$(which gatk3)
GATK3=${FileName%/*/*}/opt/gatk-3.8/GenomeAnalysisTK.jar

echo "#@##################################################################"
echo "#@ BASE QUALITY SCORE RECALIBRATION: SAMPLE ${i}"
echo
D1=`date "+%D    %T"`
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

java -Xmx80G -Djava.io.tmpdir=$TMP -jar $GATK3 \
   -T BaseRecalibrator -nct 8 \
   -R $REFERENCE \
   -I $ALIGN_MDUP/${i}_trim_"$ALIGNER"_sorted_RG_mdup.bam \
   -knownSites:vcf $VARIANTS \
   --bqsrBAQGapOpenPenalty 45 \
   -o $ALIGN_BQSR/${i}_recalibration.table

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

java -Xmx80G -Djava.io.tmpdir=$TMP -jar $GATK3 \
   -T PrintReads -nct 8 \
   -R $REFERENCE \
   -I $ALIGN_MDUP/${i}_trim_"$ALIGNER"_sorted_RG_mdup.bam \
   -BQSR $ALIGN_BQSR/${i}_recalibration.table \
   -o $ALIGN_BQSR/${i}_trim_"$ALIGNER"_sorted_RG_mdup_bqsr.bam
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
