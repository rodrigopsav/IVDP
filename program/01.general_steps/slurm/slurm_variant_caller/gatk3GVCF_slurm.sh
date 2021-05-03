#!/bin/bash

#eval "$(conda shell.bash hook)"
#conda activate ivdp
#conda activate --stack gatk3

FileName=$(which gatk3)
GATK3=${FileName%/*/*}/opt/gatk-3.8/GenomeAnalysisTK.jar

export VCF_TMP=$(echo $VCF_RAW)/"$ALIGNER"_"$CALLER"_"${ANALYSIS_ID}"
mkdir -p $VCF_TMP > /dev/null 2>&1

echo "#@##################################################################"
echo "#@ GATK 3.8 HAPLOTYPECALLER: GVCF MODE"
echo
D1=`date "+%D    %T"`
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

#https://pmbio.org/module-04-germline/0004/02/01/Germline_SnvIndel_Calling/
#https://gatkforums.broadinstitute.org/gatk/discussion/8690/combinegvcfs-input-multiple-files

java -Xmx8G -Djava.io.tmpdir=$TMP -jar $GATK3 \
-T HaplotypeCaller -nct 8 \
-R $REFERENCE \
-I $ALIGN_BQSR/"${i}"_trim_"$ALIGNER"_sorted_RG_mdup_bqsr.bam \
-o $VCF_TMP/"${i}"_"${chr}".g.vcf \
-L $chr \
-ERC GVCF \
--sample_ploidy 2 \
-variant_index_type LINEAR \
-variant_index_parameter 128000
wait

echo

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ GATK 3.8 HAPLOTYPECALLER (GVCF MODE) TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "@#############################################################"
echo




