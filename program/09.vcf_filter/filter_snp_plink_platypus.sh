#!/bin/bash 


echo "@#############################################################"
echo "@ FILTERING VCF: SNP"
echo
D1=$(date "+%D    %T")
echo "@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

echo "##### VCFTOOLS SNP HARD-FILTERING AND SELECT PASS VARIANTS"

plink2 \
--vcf $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz \
--chr-set $CHROM_SET no-xy no-mt \
--snps-only \
--min-alleles 2 \
--max-alleles 2 \
--vcf-min-qual 30 \
--vcf-min-dp $MIN_DEPTH --vcf-max-dp 250 \
--set-missing-var-ids @:# \
--maf $MAF \
--geno $MISSING \
--keep-allele-order \
--const-fid 0 \
--recode vcf id-paste=iid \
--out $VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp
wait
#--id-delim " " \
#--mind $MISSING \


# Grep header of former vcf file to avoid problems with bcftools isec and vcf-compare
zcat $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz | grep "^#" > $VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp.header
grep -v "^#" $VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp.vcf > $VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp.tmp
cat $VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp.header $VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp.tmp > $VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp.vcf
wait

rm -r $VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp.log 2> /dev/null
rm -r $VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp.header 2> /dev/null
rm -r $VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp.tmp 2> /dev/null
rm -r $VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp.nosex 2> /dev/null
wait

bgzip -f -@ $THREADS $VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp.vcf
bcftools index -f -t --threads $THREADS $VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp.vcf.gz
wait

echo
echo "##### STATS SNP PLINK2 SUMMARY #####"


### MINOR ALLELE FREQUENCY (MAF)
plink2 --vcf $VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp.vcf.gz \
--chr-set $CHROM_SET no-xy no-mt \
--snps-only \
--min-alleles 2 \
--max-alleles 2 \
--set-missing-var-ids @:# \
--freq \
--keep-allele-order \
--const-fid 0 \
--out $STATS_VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp
wait
#--freq alt1bins=0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5 \
rm -r $STATS_VCF_FILTERED_SNP/*.log 2> /dev/null
rm -r $STATS_VCF_FILTERED_SNP/*.nosex 2> /dev/null
wait
  
### MISSING DATA PER INDIVIDUAL
plink2 --vcf $VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp.vcf.gz \
--chr-set $CHROM_SET no-xy no-mt \
--snps-only \
--min-alleles 2 \
--max-alleles 2 \
--set-missing-var-ids @:# \
--missing \
--keep-allele-order \
--const-fid 0 \
--out $STATS_VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp
wait
rm -r $STATS_VCF_FILTERED_SNP/*.log 2> /dev/null
rm -r $STATS_VCF_FILTERED_SNP/*.nosex 2> /dev/null
wait

### INBREEDING
plink --vcf $VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp.vcf.gz \
--chr-set $CHROM_SET no-xy no-mt \
--snps-only --biallelic-only \
--set-missing-var-ids @:# \
--het \
--keep-allele-order \
--const-fid 0 \
--out $STATS_VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp
wait
rm -r $STATS_VCF_FILTERED_SNP/*.log 2> /dev/null
rm -r $STATS_VCF_FILTERED_SNP/*.nosex 2> /dev/null
wait

### BCFTOOLS STAT
bcftools stats $VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp.vcf.gz > $STATS_VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER".stats
wait

#plot-vcfstats -P -p $STATS_VCF_FILTERED_SNP/plot/"$ALIGNER"_"$CALLER" -t "$ALIGNER"_"$CALLER" $STATS_VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER".stats 2>/dev/null
#mv $STATS_VCF_FILTERED_SNP/plot/"$ALIGNER"_"$CALLER"/depth.0.png $STATS_VCF_FILTERED_SNP/plot/"$ALIGNER"_"$CALLER"_depth.png 2>/dev/null
#mv $STATS_VCF_FILTERED_SNP/plot/"$ALIGNER"_"$CALLER"/substitutions.0.png $STATS_VCF_FILTERED_SNP/plot/"$ALIGNER"_"$CALLER"_substitutions.png 2>/dev/null
#mv $STATS_VCF_FILTERED_SNP/plot/"$ALIGNER"_"$CALLER"/tstv_by_qual.0.png $STATS_VCF_FILTERED_SNP/plot/"$ALIGNER"_"$CALLER"_tstv_by_qual.png 2>/dev/null
#rm -r $STATS_VCF_FILTERED_SNP/plot/"$ALIGNER"_"$CALLER" 2>/dev/null
#wait

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "@ FILTERING SNP TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "@#############################################################"

echo



