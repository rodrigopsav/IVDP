#!/bin/bash

echo
echo
echo "#@#############################################################"
echo "#@               VARIANT CALLING: $OUTPUT_NAME "
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo

START1=$(date +%s)

#######################################
##### VARIANT CALLING FILES: VCFs #####
#######################################

##### Check Files from the Previous Step
for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
   source "$IVDP_DIR"/program/01.general_steps/folders.sh

   if [[ "$VARIANTS" == "none" ]]; then
      if [[ "$(ls $ALIGN_MDUP/*mdup.bam 2>/dev/null | wc -l)" -eq "0" ]]; then
         echo "#@ ERROR: No mdup.bam files from mark duplicated reads step in $ALIGN_MDUP "
         echo "Check if STEP_MDUP=yes and the log mdup.txt files in $LOG_ALIGN_MDUP"
         echo "Aborting analysis"
         exit 1      
      else
      
         for sample in $(ls $ALIGN_MDUP/*mdup.bam); do
            #if [[ "$(du $sample | awk '{print $1}')" -lt "1000" ]]; then
            if [[ "$(samtools view $sample | head -n 1000 | wc -l)" == 0 ]]; then
               echo "#@ ERROR: Malformed mdup.bam files from mark duplicated reads step in $ALIGN_MDUP: $sample"
               echo "Check if STEP_MDUP=yes and the log mdup.txt files in $LOG_ALIGN_MDUP"
               echo "Aborting analysis"
               exit 1
            fi
         done
      fi

   else

      if [[ "$(ls $ALIGN_BQSR/*bqsr.bam 2>/dev/null | wc -l)" -eq "0" ]]; then
         echo "#@ ERROR: No bqsr.bam files from base quality score recalibration step in $ALIGN_BQSR "
         echo "Check if STEP_BQSR=yes and the bqsr.txt files in $LOG_ALIGN_BQSR"
         echo "Aborting analysis"
         exit 1
      else
    
         for sample in $(ls $ALIGN_BQSR/*bqsr.bam); do
            #if [[ "$(du $sample | awk '{print $1}')" -lt "1000" ]]; then
            if [[ "$(samtools view $sample | head -n 1000 | wc -l)" == 0 ]]; then 
               echo "#@ ERROR: Malformed bqsr.bam files from base quality score recalibration step in $ALIGN_BQSR: $sample "
               echo "Check if STEP_BQSR=yes and the log bqsr.txt files in $LOG_ALIGN_BQSR"
               echo "Aborting analysis"
               exit 1
            fi
         done
      fi
      
   fi
done
wait


##### Check Output Folders for This Step
for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
   source "$IVDP_DIR"/program/01.general_steps/folders.sh
  
   mkdir -p "$VCF_RAW" > /dev/null 2>&1
   mkdir -p "$LOG_VCF_RAW" > /dev/null 2>&1
   mkdir -p "$STATS_VCF_RAW" > /dev/null 2>&1
   mkdir -p "$REPORT_VCF_RAW" > /dev/null 2>&1

done
wait


##### Run Variant Calling

### CALL_BY_CHROM == no
if [[ "${CALL_BY_CHROM,,}" == 'no' ]]; then

   n=0
   for ALIGNER in $ALIGNER_PROGRAM; do
      export ALIGNER=$ALIGNER
      source "$IVDP_DIR"/program/01.general_steps/folders.sh
      echo "##### $ALIGNER - $CALLER #####" 

      $IVDP_DIR/program/07.variant_caller/"$CALLER".sh > $LOG_VCF_RAW/"$ALIGNER"_"$CALLER".txt 2>&1 &&
      wait
         
      # Add rsIDs to VCF
      if [[ "$VARIANTS" != "none" ]]; then
         D=$(date "+%D    %T")
         echo "# Add rsIDs: $D #"
         bcftools annotate --threads $THREADS -a $VARIANTS -c ID -Oz -o $VCF_RAW/"$ALIGNER"_"$CALLER"_raw_tmp.vcf.gz $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz 
         mv $VCF_RAW/"$ALIGNER"_"$CALLER"_raw_tmp.vcf.gz $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz
         bcftools index -f -t --threads $THREADS $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz
      fi &&
      wait

      # Add chrom:pos to missing IDs
      D=$(date "+%D    %T")
      echo "# Add chrom:pos to missing IDs: $D #"
      bcftools annotate --threads $THREADS --set-id +'%CHROM\:%POS' \
      -Oz -o $VCF_RAW/"$ALIGNER"_"$CALLER"_raw_tmp.vcf.gz \
      $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz
      mv $VCF_RAW/"$ALIGNER"_"$CALLER"_raw_tmp.vcf.gz $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz
      bcftools index -f -t --threads $THREADS $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz 2> /dev/null &&
      wait

      # Annotation with snpEff
      if [[ "$ANNOTATION" != "none" ]]; then
         D=$(date "+%D    %T")
         echo "# Annotation with snpEff: $D #"
         $IVDP_DIR/program/08.vcf_annotation/snpEff_vcf.sh > /dev/null 2>&1
      fi &
      wait

      # limit jobs
      if (( $(($((++n)) % $BATCH)) == 0 )) ; then
      wait # wait until all have finished (not optimal, but most times good enough)
      echo $n Files completed
      fi
      
   done
   wait
   
   
### CALL_BY_CHROM == yes
else

   for ALIGNER in $ALIGNER_PROGRAM; do
      export ALIGNER=$ALIGNER
      source "$IVDP_DIR"/program/01.general_steps/folders.sh
      echo "##### $ALIGNER - $CALLER #####" 
      
      $IVDP_DIR/program/07.variant_caller/by_chrom/"$CALLER"_bychrom.sh > $LOG_VCF_RAW/"$ALIGNER"_"$CALLER".txt 2>&1
      wait
      
      # Add rsIDs to VCF 
      if [[ "$VARIANTS" != "none" ]]; then
         D=$(date "+%D    %T")
         echo "# Add rsIDs: $D #"
         bcftools annotate --threads $THREADS -a $VARIANTS -c ID -Oz -o $VCF_RAW/"$ALIGNER"_"$CALLER"_raw_tmp.vcf.gz $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz 
         mv $VCF_RAW/"$ALIGNER"_"$CALLER"_raw_tmp.vcf.gz $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz
         bcftools index -f -t --threads $THREADS $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz
      fi
      wait

      # Add chrom:pos to missing IDs
      D=$(date "+%D    %T")
      echo "# Add chrom:pos to missing IDs: $D #"
      bcftools annotate --threads $THREADS --set-id +'%CHROM\:%POS' \
      -Oz -o $VCF_RAW/"$ALIGNER"_"$CALLER"_raw_tmp.vcf.gz \
      $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz
      mv $VCF_RAW/"$ALIGNER"_"$CALLER"_raw_tmp.vcf.gz $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz
      bcftools index -f -t --threads $THREADS $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz 2> /dev/null
      wait
            
      # Annotation with snpEff
      if [[ "$ANNOTATION" != "none" ]]; then
         D=$(date "+%D    %T")
         echo "# Annotation with snpEff: $D #"
         $IVDP_DIR/program/08.vcf_annotation/snpEff_vcf.sh > /dev/null 2>&1
      fi
      wait
   
   done
   wait
fi
wait


###################################
######## STATISTICS RAW VCF #######
###################################
echo
D=$(date "+%D    %T")
echo "STATISTICS OF RAW VCF FILES: $D"

for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
   echo "##### $ALIGNER - $CALLER #####"
   source "$IVDP_DIR"/program/01.general_steps/folders.sh  
   
   plink2 --vcf $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz \
   --chr-set $CHROM_SET no-xy no-mt \
   --snps-only \
   --min-alleles 2 \
   --max-alleles 2 \
   --set-missing-var-ids @:# \
   --freq --missing --keep-allele-order --const-fid 0 \
   --out $STATS_VCF_RAW/"$ALIGNER"_"$CALLER" > /dev/null 2>&1 
   wait
   
   plink --vcf $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz \
   --chr-set $CHROM_SET no-xy no-mt \
   --snps-only --biallelic-only \
   --set-missing-var-ids @:# \
   --het --keep-allele-order --const-fid 0 \
   --out $STATS_VCF_RAW/"$ALIGNER"_"$CALLER" > /dev/null 2>&1 
   wait
   
   rm -r $STATS_VCF_RAW/*.log 2> /dev/null
   rm -r $STATS_VCF_RAW/*.nosex 2> /dev/null
   wait

   for chrom in $(cat $CHROM_NAME); do 
      bcftools stats --threads $THREADS -r $chrom $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz > $TMP/rawVCF_"$ALIGNER"_"$CALLER"_chrom"$chrom".stats & 
   done
   wait
   rm $TMP/rawVCF_"$ALIGNER"_"$CALLER".stats 2> /dev/null
   for chrom in $(cat $CHROM_NAME); do
      echo -e $chrom "\t" snp "\t" $(grep -w "number of SNPs:" $TMP/rawVCF_"$ALIGNER"_"$CALLER"_chrom"$chrom".stats | awk -F":" '{print $2}' | sed "s/[[:space:]]\+//g") >> $TMP/rawVCF_"$ALIGNER"_"$CALLER".stats
      echo -e $chrom "\t" indel "\t" $(grep -w "number of indels:" $TMP/rawVCF_"$ALIGNER"_"$CALLER"_chrom"$chrom".stats | awk -F":" '{print $2}' | sed "s/[[:space:]]\+//g") >> $TMP/rawVCF_"$ALIGNER"_"$CALLER".stats
   done
   wait
   unset chrom

   Rscript -e "inputPath=\"${STATS_VCF_RAW}\";aligner=\"${ALIGNER}\";caller=\"${CALLER}\";step=\"rawVCF\";chromNames=\"${CHROM_NAME}\";file1=\"${TMP}/rawVCF_${ALIGNER}_${CALLER}.stats\";file2=\"$STATS_VCF_RAW/${ALIGNER}_${CALLER}.afreq\";file3=\"$STATS_VCF_RAW/${ALIGNER}_${CALLER}.vmiss\";file4=\"$STATS_VCF_RAW/${ALIGNER}_${CALLER}.smiss\";file5=\"$STATS_VCF_RAW/${ALIGNER}_${CALLER}.het\";source(\"$IVDP_DIR/program/01.general_steps/plot_vcf_stats.R\")"
   wait
   bcftools query -f '%CHROM\t%POS\t%CHROM:%POS\n' $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz > $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.pos
   wait
   Rscript -e "inputFile=\"$VCF_RAW/${ALIGNER}_${CALLER}_raw.pos\";aligner=\"${ALIGNER}\";caller=\"${CALLER}\";outputPath=\"$STATS_VCF_RAW\";source(\"$IVDP_DIR/program/01.general_steps/plot_snpDensity.R\")"
   wait
   mv $STATS_VCF_RAW/Col1.Col0.SNP-Density.plot_DensitySNPs_"$ALIGNER"_"$CALLER".pdf $STATS_VCF_RAW/plot_DensitySNPs_"$ALIGNER"_"$CALLER"_raw.pdf
   wait

   gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dAutoRotatePages=/None -sOutputFile=$REPORT_VCF_RAW/plot_stats_"$ALIGNER"_"$CALLER"_raw.pdf $STATS_VCF_RAW/plot_stats_"$ALIGNER"_"$CALLER"_page1.pdf $STATS_VCF_RAW/plot_DensitySNPs_"$ALIGNER"_"$CALLER"_raw.pdf $STATS_VCF_RAW/plot_stats_"$ALIGNER"_"$CALLER"_page2.pdf $STATS_VCF_RAW/plot_stats_"$ALIGNER"_"$CALLER"_page3.pdf
   wait
   rm $STATS_VCF_RAW/plot_DensitySNPs_"$ALIGNER"_"$CALLER"_raw.pdf $STATS_VCF_RAW/plot_stats_"$ALIGNER"_"$CALLER"_page* $STATS_VCF_RAW/Rplots.pdf $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.pos 2> /dev/null
   wait

   bcftools stats $VCF_RAW/"$ALIGNER"_"$CALLER"_raw.vcf.gz > $STATS_VCF_RAW/"$ALIGNER"_"$CALLER".stats
   
   #plot-vcfstats -P -p $STATS_VCF_RAW/plot/"$ALIGNER"_"$CALLER" -t "$ALIGNER"_"$CALLER" $STATS_VCF_RAW/"$ALIGNER"_"$CALLER".stats 2>/dev/null
   #mv $STATS_VCF_RAW/plot/"$ALIGNER"_"$CALLER"/depth.0.png $STATS_VCF_RAW/plot/"$ALIGNER"_"$CALLER"_depth.png 2>/dev/null
   #mv $STATS_VCF_RAW/plot/"$ALIGNER"_"$CALLER"/indels.0.png $STATS_VCF_RAW/plot/"$ALIGNER"_"$CALLER"_indels.png 2>/dev/null
   #mv $STATS_VCF_RAW/plot/"$ALIGNER"_"$CALLER"/substitutions.0.png $STATS_VCF_RAW/plot/"$ALIGNER"_"$CALLER"_substitutions.png 2>/dev/null
   #mv $STATS_VCF_RAW/plot/"$ALIGNER"_"$CALLER"/tstv_by_qual.0.png $STATS_VCF_RAW/plot/"$ALIGNER"_"$CALLER"_tstv_by_qual.png 2>/dev/null
   #rm -r $STATS_VCF_RAW/plot/"$ALIGNER"_"$CALLER" 2>/dev/null
   #wait
     
done
wait

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "@ VARIANT CALLING (VCF) TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "@#############################################################"
