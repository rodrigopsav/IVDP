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
   rm -r bamlist_vs.list 2> /dev/null
   ls *bqsr.bam 2> /dev/null > bamlist_vs.list

elif [[ "$mdup" != 0 ]] && [[ "$bqsr" == 0 ]]; then
   echo "NUMBER OF MDUP BAM FILES: $mdup"
   BAM_TMP=$ALIGN_MDUP
   cd $BAM_TMP
   rm -r bamlist_vs.list 2> /dev/null
   ls *RG_mdup.bam 2> /dev/null > bamlist_vs.list
   
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
for sample in $(cat bamlist_vs.list); do
   samtools view -h $sample | head -10000 | \
   grep "@RG" | awk -F"\t" '{ print $5 }' | \
   awk -F":" '{ print $2 }' >> samples_name_vs.list
done

for chr in $(cat $CHROM_NAME); do
   samtools mpileup -B -f $REFERENCE -b bamlist_vs.list -r $chr | \
   varscan mpileup2snp \
   --min-avg-qual 20 \
   --p-value 0.05 \
   --output-vcf 1 \
   --vcf-sample-list samples_name_vs.list > $VCF_TMP/"$ALIGNER"_"$CALLER"_${chr}_raw.vcf &

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

# cd $BAM_TMP
# for chr in $(cat $CHROM_NAME); do
#    samtools mpileup -B -f $REFERENCE -b bamlist_vs.list | \
#    varscan mpileup2indel \
#    --min-avg-qual 20 \
#    --p-value 0.10 \
#    --output-vcf 1 \
#    --vcf-sample-list samples_name_vs.list > $VCF_TMP/"$ALIGNER"_"$CALLER"_${chr}_raw_indel.vcf &

   	 # limit jobs
#	 if (( $(($((++n)) % $((BATCH + 30)) )) == 0 )) ; then
#	 wait # wait until all have finished (not optimal, but most times good enough)
#	 echo $n Files completed
#	 fi
	  
# done
# wait


# Get head of vcf
# cat $(ls $VCF_TMP/*indel.vcf | head -1) | grep "^#" > $VCF_TMP/"$ALIGNER"_"$CALLER"_raw_indel.vcf
#wait

# Join vcfs
# for chr in $(cat $CHROM_NAME); do
#    grep -v '^#' $VCF_TMP/"$ALIGNER"_"$CALLER"_${chr}_raw_indel.vcf >> $VCF_TMP/"$ALIGNER"_"$CALLER"_raw_indel.vcf
# done
# wait

# bgzip -f -@ $THREADS $VCF_TMP/"$ALIGNER"_"$CALLER"_raw_indel.vcf
# bcftools index -f -t --threads $THREADS $VCF_TMP/"$ALIGNER"_"$CALLER"_raw_indel.vcf.gz
# wait

mv $VCF_TMP/"$ALIGNER"_"$CALLER"_raw.vcf* $VCF_RAW
rm -r $VCF_TMP 2> /dev/null
rm $BAM_TMP/bamlist_vs.list 2> /dev/null
rm $BAM_TMP/samples_name_vs.list 2> /dev/null

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "@ VARSCAN TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "@##################################################################"

echo
echo
echo


