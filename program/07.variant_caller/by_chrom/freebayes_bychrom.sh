#!/bin/bash

#eval "$(conda shell.bash hook)"
#conda activate ivdp
#conda activate --stack ivdp2

export VCF_TMP=$(echo $VCF_RAW)/"$ALIGNER"_"$CALLER"_"${ANALYSIS_ID}"
mkdir -p $VCF_TMP > /dev/null 2>&1

echo "#@##################################################################"
echo "#@ FREEBAYES CALLER"
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
   rm -r bamlist_fb.list 2> /dev/null
   ls *RG_mdup.bam 2> /dev/null > bamlist_fb.list
   
else
   echo "#@ ERROR: THERE ARE NO MDUP BAM FILES"
   echo "#@ PLEASE CHECK IF THE BAM FILES FROM ALIGNMENT STEP ARE OK (*RG.bam files in $ALIGN SUBFOLDERS)"
   echo "#@ PLEASE CHECK THE LOG FILES IN:"
   echo "$LOG_ALIGN"
   echo "$LOG_ALIGN_MDUP"
   exit 1
fi
echo "@##################################################################"

cd $BAM_TMP
for chr in $(cat $CHROM_NAME); do
   freebayes -L bamlist_fb.list -f $REFERENCE \
   --min-base-quality 20 --genotype-qualities \
   -r $chr -v $VCF_TMP/"$ALIGNER"_"$CALLER"_${chr}_raw.vcf &

   # limit jobs
   if (( $(($((++n)) % $((BATCH + 30)) )) == 0 )) ; then
   wait # wait until all have finished (not optimal, but most times good enough)
   echo $n Files completed
   fi

done
wait

# Get head of vcf
cat $(ls $VCF_TMP/*.vcf | head -1) | grep "^#" > $VCF_TMP/"$ALIGNER"_"$CALLER"_raw.vcf
wait

# Join vcfs
for chr in $(cat $CHROM_NAME); do
   grep -v '^#' $VCF_TMP/"$ALIGNER"_"$CALLER"_${chr}_raw.vcf >> $VCF_TMP/"$ALIGNER"_"$CALLER"_raw.vcf
done
wait

bgzip -f -@ $THREADS $VCF_TMP/"$ALIGNER"_"$CALLER"_raw.vcf
bcftools index -f -t --threads $THREADS $VCF_TMP/"$ALIGNER"_"$CALLER"_raw.vcf.gz
wait

mv $VCF_TMP/"$ALIGNER"_"$CALLER"_raw.vcf* $VCF_RAW
rm -r $VCF_TMP 2> /dev/null
rm $BAM_TMP/bamlist_fb.list 2> /dev/null

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ FREEBAYES TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "#@##################################################################"

echo
echo
echo
