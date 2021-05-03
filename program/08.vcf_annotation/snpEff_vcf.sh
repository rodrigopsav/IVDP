#!/bin/bash

### Variant annotation
#cd $VCF_RAW

snpEff -Xmx20g $ASSEMBLY -c $IVDP_DIR/dataBases/annotation_db/snpEff.config -o vcf \
   -csvStats $STATS_VCF_RAW/"$ALIGNER"_"$CALLER"_snpEff_effects.csv \
   -htmlStats $STATS_VCF_RAW/"$ALIGNER"_"$CALLER"_snpEff_effects.html \
   $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz \
   > $VCF_RAW/"$ALIGNER"_"$CALLER"_raw_ann.vcf
wait

mv $VCF_RAW/"$ALIGNER"_"$CALLER"_raw_ann.vcf $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf
bgzip -f -@ $THREADS $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf
bcftools index -f -t --threads $THREADS $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz
wait

### Feature counts from gtf file
snpEff -Xmx20g count $ASSEMBLY -c $IVDP_DIR/dataBases/annotation_db/snpEff.config \
   -noLog \
   -n $STATS_VCF_RAW/"$ALIGNER"_"$CALLER"_snpeff_counts \
   $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz
wait

### Add rsIDs
# SnpSift annotate -id $VARIANTS \
#   $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz \
#   > $VCF_RAW/"$ALIGNER"_"$CALLER"_raw_ann.vcf
# wait
# mv $VCF_RAW/"$ALIGNER"_"$CALLER"_raw_ann.vcf $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf
# bgzip -f -@ $THREADS $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf
# bcftools index -f -t --threads $THREADS $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz
#wait

multiqc -f -n "$ALIGNER"_"$CALLER"_snpEff_report.html $STATS_VCF_RAW/"$ALIGNER"_"$CALLER"_snpEff_summary.csv -o $REPORT_VCF_RAW/  > /dev/null 2>&1
wait


