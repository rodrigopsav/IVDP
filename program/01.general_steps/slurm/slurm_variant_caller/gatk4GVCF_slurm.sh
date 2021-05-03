#!/bin/bash

export VCF_TMP=$(echo $VCF_RAW)/"$ALIGNER"_"$CALLER"_"${ANALYSIS_ID}"
mkdir -p $VCF_TMP > /dev/null 2>&1

echo "#@##################################################################"
echo "#@ GATK HAPLOTYPECALLER: GVCF MODE"
echo
D1=`date "+%D    %T"`
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

#https://pmbio.org/module-04-germline/0004/02/01/Germline_SnvIndel_Calling/
#https://gatkforums.broadinstitute.org/gatk/discussion/8690/combinegvcfs-input-multiple-files

gatk --java-options '-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=4' HaplotypeCaller \
-R $REFERENCE \
-I $ALIGN_BQSR/"${i}"_trim_"$ALIGNER"_sorted_RG_mdup_bqsr.bam  \
-pairHMM LOGLESS_CACHING \
--native-pair-hmm-threads $THREADS \
-O $VCF_TMP/"${i}"_"${chr}".g.vcf \
-L $chr \
-ERC GVCF \
--sample-ploidy 2 \
--tmp-dir $TMP
wait

echo

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ GATK HAPLOTYPECALLER (GVCF MODE) TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "@#############################################################"
echo
