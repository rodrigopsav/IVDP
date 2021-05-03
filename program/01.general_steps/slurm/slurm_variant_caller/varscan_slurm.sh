#!/bin/bash

export VCF_TMP=$(echo $VCF_RAW)/"$ALIGNER"_"$CALLER"_"${ANALYSIS_ID}"
mkdir -p $VCF_TMP > /dev/null 2>&1

echo "@##################################################################"
echo "@ VARSCAN CALLER"
echo
D1=`date "+%D    %T"`
echo "@ Date and Time: $D1"
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
   rm -r bamlist_vs_"$chr".list 2> /dev/null
   ls *bqsr.bam 2> /dev/null > bamlist_vs_"$chr".list

elif [[ "$mdup" != 0 ]] && [[ "$bqsr" == 0 ]]; then
   echo "NUMBER OF MDUP BAM FILES: $mdup"
   BAM_TMP=$ALIGN_MDUP
   cd $BAM_TMP
   rm -r bamlist_vs_"$chr".list 2> /dev/null
   ls *RG_mdup.bam 2> /dev/null > bamlist_vs_"$chr".list
   
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

cd $BAM_TMP
for sample in $(cat bamlist_vs_"$chr".list); do
   samtools view -h $sample | head -10000 | \
   grep "@RG" | awk -F"\t" '{ print $5 }' | \
   awk -F":" '{ print $2 }' >> samples_name_vs_"$chr".list
done
wait

samtools mpileup -B -f $REFERENCE -b bamlist_vs_"$chr".list -r $chr | \
varscan mpileup2snp \
--min-avg-qual 20 \
--p-value 0.05 \
--output-vcf 1 \
--vcf-sample-list samples_name_vs.list > $VCF_TMP/"$ALIGNER"_"$CALLER"_${chr}_raw.vcf
wait


# cd $BAM_TMP
#samtools mpileup -B -f $REFERENCE -b bamlist_vs_"$chr".list | \
#varscan mpileup2indel \
#--min-avg-qual 20 \
#--p-value 0.10 \
#--output-vcf 1 \
#--vcf-sample-list samples_name_vs.list > $VCF_TMP/"$ALIGNER"_"$CALLER"_${chr}_raw_indel.vcf
#wait

rm $BAM_TMP/bamlist_vs_"$chr".list
rm $BAM_TMP/samples_name_vs_"$chr".list
wait

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "@ VARSCAN CHR$chr TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "@##################################################################"

echo
echo
echo


