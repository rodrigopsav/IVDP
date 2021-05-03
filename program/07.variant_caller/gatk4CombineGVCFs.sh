#!/bin/bash


echo "#@##################################################################"
echo "#@ GATK HAPLOTYPECALLER: CombineGVCFs"
echo
D1=`date "+%D    %T"`
echo "#@ Date and Time: $D1"
echo
       #################################

echo "@##### CREATE LIST WITH BAM FILES #####"
bqsr="$(ls $ALIGN_BQSR/*bqsr.bam 2> /dev/null | wc -l)"

if [[ "$bqsr" != 0 ]]; then
   echo "NUMBER OF BQSR BAM FILES: $bqsr"
   cd $ALIGN_BQSR
   rm -r bamlist_GATK.list 2> /dev/null
   ls *bqsr.bam 2> /dev/null > bamlist_GATK.list
   
else
   echo "#@ ERROR: THERE ARE NO BQSR BAM FILES"
   echo "#@ PLEASE CHECK IF THE BAM FILES FROM ALIGNMENT STEP ARE OK (*RG.bam files in $ALIGN SUBFOLDERS)"
   echo "#@ PLEASE CHECK THE LOG FILES IN:"
   echo "$LOG_ALIGN"
   echo "$LOG_ALIGN_MDUP"
   echo "$LOG_ALIGN_BQSR"
   exit 1
fi
echo "@##################################################################"

echo
echo "#####@ GATK CombineGVCFs #####"
START1=$(date +%s)

#https://pmbio.org/module-04-germline/0004/02/01/Germline_SnvIndel_Calling/
#https://gatkforums.broadinstitute.org/gatk/discussion/8690/combinegvcfs-input-multiple-files

cd $GVCF_RAW/"$ALIGNER"_"$CALLER"
ls *g.vcf > input_gvcf.list
wait

n=0
for chr in $(cat $CHROM_NAME); do
   gatk --java-options '-Xmx16g' CombineGVCFs \
   -R $REFERENCE \
   -V $GVCF_RAW/"$ALIGNER"_"$CALLER"/input_gvcf.list \
   -O $GVCF_RAW/"$ALIGNER"_"$CALLER"/"$ALIGNER"_"$CALLER"_"$chr"_combined.g.vcf \
   -L $chr \
   --tmp-dir $TMP &

	# limit jobs
	if (( $(($((++n)) % $((BATCH + 30)) )) == 0 )) ; then
	wait # wait until all have finished (not optimal, but most times good enough)
	echo $n Files completed
	fi
done
wait

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#####@ GATK GenotypeGVCFs #####"
START2=$(date +%s)

#https://pmbio.org/module-04-germline/0004/02/01/Germline_SnvIndel_Calling/

cd $GVCF_RAW/"$ALIGNER"_"$CALLER"
n=0
for chr in $(cat $CHROM_NAME); do
   gatk --java-options '-Xmx16g -XX:+UseParallelGC -XX:ParallelGCThreads=4' GenotypeGVCFs \
   -R $REFERENCE \
   -V $GVCF_RAW/"$ALIGNER"_"$CALLER"/"$ALIGNER"_"$CALLER"_"$chr"_combined.g.vcf \
   -O $GVCF_RAW/"$ALIGNER"_"$CALLER"/"$ALIGNER"_"$CALLER"_"$chr"_raw.vcf \
   -L $chr \
   --sample-ploidy 2 \
   --tmp-dir $TMP &

	# limit jobs
	if (( $(($((++n)) % $((BATCH + 30)) )) == 0 )) ; then
	wait # wait until all have finished (not optimal, but most times good enough)
	echo $n Files completed
	fi
done
wait

# Get head of vcf
cat $(ls $GVCF_RAW/"$ALIGNER"_"$CALLER"/*raw.vcf | head -1) | grep "^#" > $GVCF_RAW/"$ALIGNER"_"$CALLER"/"$ALIGNER"_"$CALLER"_raw.vcf
wait

# Join vcfs
for chr in $(cat $CHROM_NAME); do
   grep -v '^#' $GVCF_RAW/"$ALIGNER"_"$CALLER"/"$ALIGNER"_"$CALLER"_"${chr}"_raw.vcf >> $GVCF_RAW/"$ALIGNER"_"$CALLER"/"$ALIGNER"_"$CALLER"_raw.vcf
done
wait

bgzip -f -@ $THREADS $GVCF_RAW/"$ALIGNER"_"$CALLER"/"$ALIGNER"_"$CALLER"_raw.vcf
bcftools index -f -t --threads $THREADS $GVCF_RAW/"$ALIGNER"_"$CALLER"/"$ALIGNER"_"$CALLER"_raw.vcf.gz
wait

mv $GVCF_RAW/"$ALIGNER"_"$CALLER"/"$ALIGNER"_"$CALLER"_raw.vcf* $VCF_RAW
rm -r $GVCF_RAW/"$ALIGNER"_"$CALLER"/input_gvcf.list 2> /dev/null
rm -r $GVCF_RAW/"$ALIGNER"_"$CALLER"/"$ALIGNER"_"$CALLER"_*_combined.*
rm -r $GVCF_RAW/"$ALIGNER"_"$CALLER"/"$ALIGNER"_"$CALLER"_*_raw.* 2> /dev/null
rm $ALIGN_BQSR/bamlist_GATK.list 2> /dev/null

END2=$(date +%s)
DIFF2=$(( $END2 - $START2 ))
DIFF3=$(( $END2 - $START1 ))

echo
echo "#@ GATK CombineGVCFs TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@ GATK GenotypeGVCFs TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF2/3600)) $(($DIFF2%3600/60)) $(($DIFF2%60)))"
echo
echo "#@ GATK CONVERT GVCF TO VCF TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF3/3600)) $(($DIFF3%3600/60)) $(($DIFF3%60)))"
echo
echo "#@#############################################################"
echo

