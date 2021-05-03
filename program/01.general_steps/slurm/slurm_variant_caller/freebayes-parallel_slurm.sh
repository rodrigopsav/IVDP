#!/bin/bash

#eval "$(conda shell.bash hook)"
#conda activate ivdp
#conda activate --stack ivdp2

export VCF_TMP=$(echo $VCF_RAW)/"$ALIGNER"_"$CALLER"_"${ANALYSIS_ID}"
mkdir -p $VCF_TMP > /dev/null 2>&1

echo "#@##################################################################"
echo "#@ FREEBAYES-PARALLEL CALLER"
echo
D1=`date "+%D    %T"`
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

echo "@##### CREATE LIST WITH BAM FILES #####"
mdup="$(ls $ALIGN_MDUP/*RG_mdup.bam 2> /dev/null | wc -l)"

if [[ "$mdup" != 0 ]]; then
   echo "NUMBER OF MDUP BAM FILES: $mdup"
   BAM_TMP=$ALIGN_MDUP
   cd $BAM_TMP
   rm -r bamlist_fb_"$chr".list 2> /dev/null
   ls *RG_mdup.bam 2> /dev/null > bamlist_fb_"$chr".list
   
else
   echo "#@ ERROR: THERE ARE NO MDUP BAM FILES"
   echo "#@ PLEASE CHECK IF THE BAM FILES FROM ALIGNMENT STEP ARE OK (*RG.bam files in $ALIGN SUBFOLDERS)"
   echo "#@ PLEASE CHECK THE LOG FILES IN:"
   echo "$LOG_ALIGN"
   echo "$LOG_ALIGN_MDUP"
   exit 1
fi
echo "@##################################################################"


export REFERENCE_TMP=$(echo $REFERENCE | rev | cut -d '/' -f 2- | rev)
wait
samtools faidx $REFERENCE $chr > $REFERENCE_TMP/ref_"$chr".fa
wait
samtools faidx $REFERENCE_TMP/ref_"$chr".fa
wait

cd $BAM_TMP
freebayes-parallel <(fasta_generate_regions.py $REFERENCE_TMP/ref_"$chr".fa.fai 100000) \
$THREADS -f $REFERENCE_TMP/ref_"$chr".fa -L bamlist_fb_"$chr".list \
--min-base-quality 20 --genotype-qualities >$VCF_TMP/"$ALIGNER"_"$CALLER"_${chr}_raw.vcf
wait
   
rm $REFERENCE_TMP/ref_"$chr".fa 2> /dev/null
rm $REFERENCE_TMP/ref_"$chr".fa.fai 2> /dev/null
rm $BAM_TMP/bamlist_fb_"$chr".list
wait

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ FREEBAYES-PARALLEL CHR$chr TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "#@##################################################################"

echo
echo
echo


