#!/bin/bash


echo "#@##################################################################"
echo "#@ SPLIT’N’TRIM CIGAR: SAMPLE ${i}"
echo
D1=`date "+%D    %T"`
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

gatk --java-options '-Djava.io.tmpdir=$TMP -Xmx20g' SplitNCigarReads \
	-R $REFERENCE \
	-I $ALIGN_MDUP/${i}_trim_"$ALIGNER"_sorted_RG_mdup.bam \
	-O $ALIGN_MDUP/${i}_trim_"$ALIGNER"_sorted_RG_mdup_splitncigar.bam \
	--tmp-dir $TMP
wait


END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ SPLIT’N’TRIM SAMPLE $i TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "#@##################################################################"

echo

