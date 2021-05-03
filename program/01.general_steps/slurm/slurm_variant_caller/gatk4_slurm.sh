#!/bin/bash

export VCF_TMP=$(echo $VCF_RAW)/"$ALIGNER"_"$CALLER"_"${ANALYSIS_ID}"
mkdir -p $VCF_TMP > /dev/null 2>&1

echo "#@##################################################################"
echo "#@ GATK HAPLOTYPECALLER"
echo
D1=`date "+%D    %T"`
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

echo "@##### CREATE LIST WITH BAM FILES #####"
mdup="$(ls $ALIGN_MDUP/*RG_mdup.bam 2> /dev/null | wc -l)"
bqsr="$(ls $ALIGN_BQSR/*bqsr.bam 2> /dev/null | wc -l)"

if [[ "$mdup" != 0 ]] && [[ "$bqsr" != 0 ]]; then
   echo "NUMBER OF BQSR BAM FILES: $bqsr"
   BAM_TMP=$ALIGN_BQSR
   cd $BAM_TMP
   rm -r bamlist_GATK_"$chr".list 2> /dev/null
   ls *bqsr.bam 2> /dev/null > bamlist_GATK_"$chr".list

elif [[ "$mdup" != 0 ]] && [[ "$bqsr" == 0 ]]; then
   echo "NUMBER OF MDUP BAM FILES: $mdup"
   BAM_TMP=$ALIGN_MDUP
   cd $BAM_TMP
   rm -r bamlist_GATK_"$chr".list 2> /dev/null
   ls *RG_mdup.bam 2> /dev/null > bamlist_GATK_"$chr".list
   
else
   echo "#@ ERROR: THERE ARE NO MDUP OR BQSR BAM FILES"
   echo "#@ PLEASE CHECK IF THE BAM FILES FROM ALIGNMENT STEP ARE OK (*RG.bam files in $ALIGN SUBFOLDERS)"
   echo "#@ PLEASE CHECK THE LOG FILES IN:"
   echo "$LOG_ALIGN"
   echo "$LOG_ALIGN_MDUP"
   echo "$LOG_ALIGN_BQSR"
   exit 1
fi
echo "@##################################################################"

echo
echo
echo

#https://pmbio.org/module-04-germline/0004/02/01/Germline_SnvIndel_Calling/
#https://gatkforums.broadinstitute.org/gatk/discussion/8690/combinegvcfs-input-multiple-files

cd $BAM_TMP
gatk --java-options '-Xmx32g -XX:+UseParallelGC -XX:ParallelGCThreads=4' HaplotypeCaller \
-R $REFERENCE \
-I bamlist_GATK_"$chr".list \
-pairHMM LOGLESS_CACHING \
--native-pair-hmm-threads $THREADS \
-O $VCF_TMP/"$ALIGNER"_"$CALLER"_${chr}_raw.vcf \
-L $chr \
--sample-ploidy 2 \
--tmp-dir $TMP
wait

rm $BAM_TMP/bamlist_GATK_"$chr".list
wait

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ GATK HAPLOTYPECALLER CHR$chr TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "#@#############################################################"

echo
echo
echo
