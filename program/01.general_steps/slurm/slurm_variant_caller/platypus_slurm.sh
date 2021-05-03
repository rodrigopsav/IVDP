#!/bin/bash

#eval "$(conda shell.bash hook)"
#conda activate ivdp
#conda activate --stack ivdp2

export VCF_TMP=$(echo $VCF_RAW)/"$ALIGNER"_"$CALLER"_"${ANALYSIS_ID}"
mkdir -p $VCF_TMP > /dev/null 2>&1

echo "#@##################################################################"
echo "#@ PLATYPUS CALLER"
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
   rm -r bamlist_platypus_"$chr".list 2> /dev/null
   ls *bqsr.bam 2> /dev/null > bamlist_platypus_"$chr".list

elif [[ "$mdup" != 0 ]] && [[ "$bqsr" == 0 ]]; then
   echo "NUMBER OF MDUP BAM FILES: $mdup"
   BAM_TMP=$ALIGN_MDUP
   cd $BAM_TMP
   rm -r bamlist_platypus_"$chr".list 2> /dev/null
   ls *RG_mdup.bam 2> /dev/null > bamlist_platypus_"$chr".list
   
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
platypus callVariants \
--bamFiles=bamlist_platypus_"$chr".list \
--refFile=$REFERENCE \
--output=$VCF_TMP/"$ALIGNER"_"$CALLER"_${chr}_raw.vcf \
--regions="$chr" \
--nCPU=$THREADS
wait

rm $BAM_TMP/log.txt
rm $BAM_TMP/bamlist_platypus_"$chr".list
wait

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ PLATYPUS TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "#@##################################################################"

echo
echo
echo

