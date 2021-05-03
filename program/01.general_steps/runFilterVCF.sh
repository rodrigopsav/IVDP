#!/bin/bash

echo
echo
echo "#@#############################################################"
echo "#@ VCF FILTERING: $OUTPUT_NAME "
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo

START1=$(date +%s)

#########################
##### FILTERING VCF #####
#########################

##### Check Files from the Previous Step
source "$IVDP_DIR"/program/01.general_steps/folders.sh
cd $VCF_RAW
if [[ "$(ls *raw.vcf.gz 2>/dev/null | wc -l)" -eq "0" ]]; then
   echo "#@ ERROR: No raw.vcf.gz files from variant calling step in $VCF_RAW "
   echo "Check if STEP_VCF_CALL=yes and STEP_GVCF_TO_VCF=yes (if gvcf mode) and check the log files in $LOG_VCF_RAW"
   echo "Aborting analysis"
   exit 1
else
   for file in $(ls *raw.vcf.gz); do
      if [[ "$(bcftools view -H $file | head -n 1000 | wc -l )" == 0 ]]; then 
         echo "#@ ERROR: No variants in raw.vcf.gz files from variant calling step in $VCF_RAW: $file"
         echo "Check if STEP_VCF_CALL=yes and STEP_GVCF_TO_VCF=yes (if gvcf mode) and check the log files in $LOG_VCF_RAW"
         echo "Aborting analysis"
         exit 1
      fi
   done
fi
wait


##### Check Output Folders for This Step
source $IVDP_DIR/program/01.general_steps/folders.sh
mkdir -p "$VCF_FILTERED_SNP" > /dev/null 2>&1
mkdir -p "$VCF_FILTERED_INDEL" > /dev/null 2>&1
mkdir -p "$LOG_VCF_FILTERED" > /dev/null 2>&1
mkdir -p "$STATS_VCF_FILTERED_SNP" > /dev/null 2>&1
mkdir -p "$STATS_VCF_FILTERED_INDEL" > /dev/null 2>&1
mkdir -p "$REPORT_VCF_FILTERED_SNP" > /dev/null 2>&1
mkdir -p "$REPORT_VCF_FILTERED_INDEL" > /dev/null 2>&1
wait


##### Run Filter VCF
for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
   echo "##### $ALIGNER #####"

   for CALLER in $CALLER_PROGRAM; do
      export CALLER=$CALLER
      echo "$CALLER"
      source "$IVDP_DIR"/program/01.general_steps/folders.sh

      #if [[ $CALLER == "platypus" ]]; then
      #   $IVDP_DIR/program/09.filter_vcf/filter_snp_vcftools_platypus.sh > $LOG_VCF_FILTERED/"$ALIGNER"_"$CALLER"_filter_snp.txt 2>&1 &
      #   $IVDP_DIR/program/09.filter_vcf/filter_indel_vcftools_platypus.sh > $LOG_VCF_FILTERED/"$ALIGNER"_"$CALLER"_filter_indel.txt 2>&1 &
      #else
      #   $IVDP_DIR/program/09.filter_vcf/filter_snp_vcftools.sh > $LOG_VCF_FILTERED/"$ALIGNER"_"$CALLER"_filter_snp.txt 2>&1 &
      #   $IVDP_DIR/program/09.filter_vcf/filter_indel_vcftools.sh > $LOG_VCF_FILTERED/"$ALIGNER"_"$CALLER"_filter_indel.txt 2>&1 &
      #fi
      #wait
      
      
      if [[ $CALLER == "platypus" ]]; then
         $IVDP_DIR/program/09.vcf_filter/filter_snp_plink_platypus.sh > $LOG_VCF_FILTERED/"$ALIGNER"_"$CALLER"_filter_snp.txt 2>&1 &
         $IVDP_DIR/program/09.vcf_filter/filter_indel_plink_platypus.sh > $LOG_VCF_FILTERED/"$ALIGNER"_"$CALLER"_filter_indel.txt 2>&1 &
      else
         $IVDP_DIR/program/09.vcf_filter/filter_snp_plink.sh > $LOG_VCF_FILTERED/"$ALIGNER"_"$CALLER"_filter_snp.txt 2>&1 &
         $IVDP_DIR/program/09.vcf_filter/filter_indel_plink.sh > $LOG_VCF_FILTERED/"$ALIGNER"_"$CALLER"_filter_indel.txt 2>&1 &
      fi
      wait

      #STATISTICS
      for chrom in $(cat $CHROM_NAME); do 
         bcftools stats --threads $THREADS -r $chrom $VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp.vcf.gz > $TMP/filteredVCF_"$ALIGNER"_"$CALLER"_chrom"$chrom".stats & 
      done
      wait
      rm $TMP/filteredVCF_"$ALIGNER"_"$CALLER".stats 2> /dev/null
      for chrom in $(cat $CHROM_NAME); do
         echo -e $chrom "\t" snp "\t" $(grep -w "number of SNPs:" $TMP/filteredVCF_"$ALIGNER"_"$CALLER"_chrom"$chrom".stats | awk -F":" '{print $2}' | sed "s/[[:space:]]\+//g") >> $TMP/filteredVCF_"$ALIGNER"_"$CALLER".stats
      done
      wait
      unset chrom
   
      Rscript -e "inputPath=\"${STATS_VCF_FILTERED_SNP}\";aligner=\"${ALIGNER}\";caller=\"${CALLER}\";step=\"filteredVCF\";chromNames=\"${CHROM_NAME}\";file1=\"${TMP}/filteredVCF_${ALIGNER}_${CALLER}.stats\";file2=\"$STATS_VCF_FILTERED_SNP/${ALIGNER}_${CALLER}_filter_snp.afreq\";file3=\"$STATS_VCF_FILTERED_SNP/${ALIGNER}_${CALLER}_filter_snp.vmiss\";file4=\"$STATS_VCF_FILTERED_SNP/${ALIGNER}_${CALLER}_filter_snp.smiss\";file5=\"$STATS_VCF_FILTERED_SNP/${ALIGNER}_${CALLER}_filter_snp.het\";source(\"$IVDP_DIR/program/01.general_steps/plot_vcf_stats.R\")"
      wait
      bcftools query -f '%CHROM\t%POS\t%CHROM:%POS\n' $VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp.vcf.gz > $VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp.pos
      wait
      Rscript -e "inputFile=\"${VCF_FILTERED_SNP}/${ALIGNER}_${CALLER}_filter_snp.pos\";aligner=\"${ALIGNER}\";caller=\"${CALLER}\";outputPath=\"${STATS_VCF_FILTERED_SNP}\";source(\"$IVDP_DIR/program/01.general_steps/plot_snpDensity.R\")"
      wait
      mv $STATS_VCF_FILTERED_SNP/Col1.Col0.SNP-Density.plot_DensitySNPs_"$ALIGNER"_"$CALLER".pdf $STATS_VCF_FILTERED_SNP/plot_DensitySNPs_"$ALIGNER"_"$CALLER"_filter_snp.pdf
      wait
      
      gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dAutoRotatePages=/None -sOutputFile=$REPORT_VCF_FILTERED_SNP/plot_stats_"$ALIGNER"_"$CALLER"_filter_snp.pdf $STATS_VCF_FILTERED_SNP/plot_stats_"$ALIGNER"_"$CALLER"_page1.pdf $STATS_VCF_FILTERED_SNP/plot_DensitySNPs_"$ALIGNER"_"$CALLER"_filter_snp.pdf $STATS_VCF_FILTERED_SNP/plot_stats_"$ALIGNER"_"$CALLER"_page2.pdf $STATS_VCF_FILTERED_SNP/plot_stats_"$ALIGNER"_"$CALLER"_page3.pdf
      wait
      rm $STATS_VCF_FILTERED_SNP/plot_DensitySNPs_"$ALIGNER"_"$CALLER"_filter_snp.pdf $STATS_VCF_FILTERED_SNP/plot_stats_"$ALIGNER"_"$CALLER"_page* $STATS_VCF_FILTERED_SNP/Rplots.pdf $VCF_FILTERED_SNP/"$ALIGNER"_"$CALLER"_filter_snp.pos 2> /dev/null
      wait
      
   done
done
wait

#################################
##### MULTIQC VCF FILTERING #####
#################################
echo
echo "MULTIQC REPORT"

multiqc -f -n report_vcf_filter_snp.html $STATS_VCF_FILTERED_SNP/* -o $REPORT_VCF_FILTERED_SNP > /dev/null 2>&1
multiqc -f -n report_vcf_filter_indel.html $STATS_VCF_FILTERED_INDEL/* -o $REPORT_VCF_FILTERED_INDEL > /dev/null 2>&1
wait


END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ VCF FILTERING TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@#############################################################"

echo
echo

echo "#@#############################################################"
echo "#@ COMBINE VCF: $OUTPUT_NAME "
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo

START2=$(date +%s)


############################
####### COMBINE VCF ########
## IF MORE THAN 2 CALLERS ##
############################


##### Check Files from the Previous Step
source "$IVDP_DIR"/program/01.general_steps/folders.sh
cd $VCF_FILTERED_SNP
if [[ "$(ls *filter_snp.vcf.gz 2>/dev/null | wc -l)" -eq "0" ]]; then
   echo "#@ ERROR: No filter_snp.vcf.gz files from filtered variant calling step in $VCF_FILTERED_SNP "
   echo "Check if STEP_VCF_FILTER=yes and the log files in $LOG_VCF_FILTERED"
   echo "Aborting analysis"
   exit 1
else
   for file in $(ls *filter_snp.vcf.gz); do
      if [[ "$(bcftools view -H $file | head -n 1000 | wc -l )" == 0 ]]; then 
         echo "#@ Warning: No variants in filter_snp.vcf.gz files from filtered variant calling step in $VCF_FILTERED_SNP: $file "
         echo "Check if STEP_VCF_FILTER=yes and the log files in $LOG_VCF_FILTERED"
         echo "Aborting analysis"
         exit 1
      fi
   done
fi
wait


##### Run Combine VCF
#https://unix.stackexchange.com/questions/499027/sum-and-count-in-for-loop
export N_CALLERS=0
for CALLER in $CALLER_PROGRAM; do
   N_CALLERS=$((N_CALLERS + 1))
done   


if [[ "$N_CALLERS" -gt "1" ]]; then

   if echo $COMBINE_VCF | grep -q "partial"; then
      echo "COMBINE VCF: PARTIAL (Consider SNPs that appeared at least in 2 different vcfs - from different variant callers nested to a aligner)"
      export SET=2
      
      ##### Create Folders for This Step
      source "$IVDP_DIR"/program/01.general_steps/folders.sh
      mkdir -p "$VCF_FILTERED_COMBINED/partial" > /dev/null 2>&1
      #mkdir -p "$LOG_VCF_FILTERED_COMBINED" > /dev/null 2>&1
      mkdir -p "$STATS_VCF_FILTERED_COMBINED/partial" > /dev/null 2>&1
      mkdir -p "$REPORT_VCF_FILTERED_COMBINED/partial" > /dev/null 2>&1
      wait

      echo
      for ALIGNER in $ALIGNER_PROGRAM; do
         export ALIGNER=$ALIGNER
         echo "##### $ALIGNER #####"
         source "$IVDP_DIR"/program/01.general_steps/folders.sh
         
         # COMPARE FILTERED SNP BY ALIGNER
         vcf-compare $VCF_FILTERED_SNP/"$ALIGNER"*filter_snp.vcf.gz \
         > $STATS_VCF_FILTERED_COMBINED/partial/"$ALIGNER"_snp.combined 2> /dev/null
         wait
         grep "^VN" $STATS_VCF_FILTERED_COMBINED/partial/"$ALIGNER"_snp.combined | cut -f 2- | \
         sed -e "s|$VCF_FILTERED_SNP||g" > $STATS_VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp.summary
         wait
         rm -r $STATS_VCF_FILTERED_COMBINED/partial/"$ALIGNER"_snp.combined 2> /dev/null
         wait
        
         # COMBINE FILTERED SNP BY ALIGNER
         bcftools isec --threads $THREADS -Oz -p $VCF_FILTERED_COMBINED/partial -n+$SET -w1 \
         $VCF_FILTERED_SNP/"$ALIGNER"*filter_snp.vcf.gz 2> /dev/null
         wait
        
         #bcftools concat --threads $THREADS -a -D \
         #-Oz -o $STATS_VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp.vcf.gz \
         #$STATS_VCF_FILTERED_COMBINED/partial/000*.vcf.gz
         #wait
        
         mv $VCF_FILTERED_COMBINED/partial/0000.vcf.gz $VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp.vcf.gz
         tabix -f -p vcf $VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp.vcf.gz
         wait
        
         rm $VCF_FILTERED_COMBINED/partial/000*
         rm $VCF_FILTERED_COMBINED/partial/README.txt
         rm $VCF_FILTERED_COMBINED/partial/sites.txt
         wait
         
      done
      wait
   fi
   wait

   #
   #
   #
   
   if echo $COMBINE_VCF | grep -q "full"; then
      echo "COMBINE VCF: FULL (Consider SNPs that appeared in all vcfs - from different variant callers nested to a aligner)"
      export SET=$N_CALLERS

      ##### Create Folders for This Step
      source "$IVDP_DIR"/program/01.general_steps/folders.sh
      mkdir -p "$VCF_FILTERED_COMBINED/full" > /dev/null 2>&1
      #mkdir -p "$LOG_VCF_FILTERED_COMBINED" > /dev/null 2>&1
      mkdir -p "$STATS_VCF_FILTERED_COMBINED/full" > /dev/null 2>&1
      mkdir -p "$REPORT_VCF_FILTERED_COMBINED/full" > /dev/null 2>&1
      wait
      
      echo
      for ALIGNER in $ALIGNER_PROGRAM; do
         export ALIGNER=$ALIGNER
         echo "##### $ALIGNER #####"
         source "$IVDP_DIR"/program/01.general_steps/folders.sh
         
         # COMPARE FILTERED SNP BY ALIGNER
         vcf-compare $VCF_FILTERED_SNP/"$ALIGNER"*filter_snp.vcf.gz \
         > $STATS_VCF_FILTERED_COMBINED/full/"$ALIGNER"_snp.combined 2> /dev/null
         wait
         grep "^VN" $STATS_VCF_FILTERED_COMBINED/full/"$ALIGNER"_snp.combined | cut -f 2- | \
         sed -e "s|$VCF_FILTERED_SNP||g" > $STATS_VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp.summary
         wait
         rm -r $STATS_VCF_FILTERED_COMBINED/full/"$ALIGNER"_snp.combined 2> /dev/null
         wait
        
         # COMBINE FILTERED SNP BY ALIGNER
         bcftools isec --threads $THREADS -Oz -p $VCF_FILTERED_COMBINED/full -n+$SET -w1 \
         $VCF_FILTERED_SNP/"$ALIGNER"*filter_snp.vcf.gz 2> /dev/null
         wait
        
         #bcftools concat --threads $THREADS -a -D \
         #-Oz -o $STATS_VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp.vcf.gz \
         #$STATS_VCF_FILTERED_COMBINED/full/000*.vcf.gz
         #wait
        
         mv $VCF_FILTERED_COMBINED/full/0000.vcf.gz $VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp.vcf.gz
         tabix -f -p vcf $VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp.vcf.gz
         wait
        
         rm $VCF_FILTERED_COMBINED/full/000*
         rm $VCF_FILTERED_COMBINED/full/README.txt
         rm $VCF_FILTERED_COMBINED/full/sites.txt
         wait
      
      done
      wait
   fi
   wait
   
fi
wait


##################################
##### STATISTICS VCF COMBINE #####
##################################
echo
echo "STATISTICS COMBINED VCF"

if [[ "$N_CALLERS" -gt "1" ]]; then

   if echo $COMBINE_VCF | grep -q "partial"; then

      for ALIGNER in $ALIGNER_PROGRAM; do
         export ALIGNER=$ALIGNER
         echo "##### $ALIGNER #####"
         source "$IVDP_DIR"/program/01.general_steps/folders.sh
              
           ### MEAN SITE DEPTH
           #vcftools --gzvcf $VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp.vcf.gz \
           #--site-mean-depth --out $STATS_VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp > /dev/null 2>&1
           #wait
           
           
           ### MINOR ALLELE FREQUENCY (MAF)
           #vcftools --gzvcf $VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp.vcf.gz \
           #--freq --out $STATS_VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp > /dev/null 2>&1
           #wait
           
           plink2 --vcf $VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp.vcf.gz \
           --chr-set $CHROM_SET no-xy no-mt \
           --snps-only \
           --min-alleles 2 \
           --max-alleles 2 \
           --set-missing-var-ids @:# \
           --freq \
           --keep-allele-order \
           --const-fid 0 \
           --out $STATS_VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp > /dev/null 2>&1
           wait
           #--freq alt1bins=0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5 \
           rm -r $STATS_VCF_FILTERED_COMBINED/partial/*.log 2> /dev/null
           rm -r $STATS_VCF_FILTERED_COMBINED/partial/*.nosex 2> /dev/null
           wait
           
           ### MISSING DATA PER INDIVIDUAL
           #vcftools --gzvcf $VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp.vcf.gz \
           #--missing-indv --out $STATS_VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp > /dev/null 2>&1
           #wait

           plink2 --vcf $VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp.vcf.gz \
           --chr-set $CHROM_SET no-xy no-mt \
           --snps-only \
           --min-alleles 2 \
           --max-alleles 2 \
           --set-missing-var-ids @:# \
           --missing \
           --keep-allele-order \
           --const-fid 0 \
           --out $STATS_VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp > /dev/null 2>&1
           wait
           rm -r $STATS_VCF_FILTERED_COMBINED/partial/*.log 2> /dev/null
           rm -r $STATS_VCF_FILTERED_COMBINED/partial/*.nosex 2> /dev/null
           wait     
           
           
           ### INBREEDING
            #vcftools --gzvcf $VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp.vcf.gz \
           #--het --out $STATS_VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp > /dev/null 2>&1
           #wait          
           
           plink --vcf $VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp.vcf.gz \
           --chr-set $CHROM_SET no-xy no-mt \
           --snps-only --biallelic-only \
           --set-missing-var-ids @:# \
           --het \
           --keep-allele-order \
           --const-fid 0 \
           --out $STATS_VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp > /dev/null 2>&1
           wait
           rm -r $STATS_VCF_FILTERED_COMBINED/partial/*.log 2> /dev/null
           rm -r $STATS_VCF_FILTERED_COMBINED/partial/*.nosex 2> /dev/null
           wait
           
           for chrom in $(cat $CHROM_NAME); do 
              bcftools stats --threads $THREADS -r $chrom $VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp.vcf.gz > $TMP/partialCombinedVCF_"$ALIGNER"_chrom"$chrom".stats & 
           done
           wait
           rm $TMP/partialCombinedVCF_"$ALIGNER".stats 2> /dev/null
           for chrom in $(cat $CHROM_NAME); do
              echo -e $chrom "\t" snp "\t" $(grep -w "number of SNPs:" $TMP/partialCombinedVCF_"$ALIGNER"_chrom"$chrom".stats | awk -F":" '{print $2}' | sed "s/[[:space:]]\+//g") >> $TMP/partialCombinedVCF_"$ALIGNER".stats
           done
           wait
           unset chrom
           
           ### STATISTICS
           Rscript -e "inputPath=\"$STATS_VCF_FILTERED_COMBINED/partial\";aligner=\"${ALIGNER}\";caller=\"none\";step=\"partialCombinedVCF\";chromNames=\"${CHROM_NAME}\";file1=\"${TMP}/partialCombinedVCF_${ALIGNER}.stats\";file2=\"$STATS_VCF_FILTERED_COMBINED/partial/${ALIGNER}_combined_snp.afreq\";file3=\"$STATS_VCF_FILTERED_COMBINED/partial/${ALIGNER}_combined_snp.vmiss\";file4=\"$STATS_VCF_FILTERED_COMBINED/partial/${ALIGNER}_combined_snp.smiss\";file5=\"$STATS_VCF_FILTERED_COMBINED/partial/${ALIGNER}_combined_snp.het\";source(\"$IVDP_DIR/program/01.general_steps/plot_vcf_stats.R\")"
           wait
           bcftools query -f '%CHROM\t%POS\t%CHROM:%POS\n' $VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp.vcf.gz > $VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp.pos
           wait
           Rscript -e "inputFile=\"$VCF_FILTERED_COMBINED/partial/${ALIGNER}_combined_snp.pos\";aligner=\"${ALIGNER}\";caller=\"none\";outputPath=\"$STATS_VCF_FILTERED_COMBINED/partial\";source(\"$IVDP_DIR/program/01.general_steps/plot_snpDensity.R\")"
           wait
           mv $STATS_VCF_FILTERED_COMBINED/partial/Col1.Col0.SNP-Density.plot_DensitySNPs_"$ALIGNER"_none.pdf $STATS_VCF_FILTERED_COMBINED/partial/plot_DensitySNPs_"$ALIGNER"_combined_snp.pdf
           wait
                     
           gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dAutoRotatePages=/None -sOutputFile=$REPORT_VCF_FILTERED_COMBINED/partial/plot_stats_"$ALIGNER"_combined_snp.pdf $STATS_VCF_FILTERED_COMBINED/partial/plot_stats_"$ALIGNER"_none_page1.pdf $STATS_VCF_FILTERED_COMBINED/partial/plot_DensitySNPs_"$ALIGNER"_combined_snp.pdf $STATS_VCF_FILTERED_COMBINED/partial/plot_stats_"$ALIGNER"_none_page2.pdf $STATS_VCF_FILTERED_COMBINED/partial/plot_stats_"$ALIGNER"_none_page3.pdf
           wait
           rm $STATS_VCF_FILTERED_COMBINED/partial/plot_DensitySNPs_"$ALIGNER"_combined_snp.pdf $STATS_VCF_FILTERED_COMBINED/partial/plot_stats_"$ALIGNER"_none_page* $STATS_VCF_FILTERED_COMBINED/partial/Rplots.pdf $VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp.pos 2> /dev/null
           wait

           bcftools stats $VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp.vcf.gz > $STATS_VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp.stats
           wait
           
           #plot-vcfstats -P -p $STATS_VCF_FILTERED_COMBINED/partial/plot/"$ALIGNER" -t ${ALIGNER}_combined_snp $STATS_VCF_FILTERED_COMBINED/partial/"$ALIGNER"_combined_snp.stats 2>/dev/null
           #mv $STATS_VCF_FILTERED_COMBINED/partial/plot/"$ALIGNER"/depth.0.png $STATS_VCF_FILTERED_COMBINED/partial/plot/"$ALIGNER"_combined_snp_depth.png 2>/dev/null
           #mv $STATS_VCF_FILTERED_COMBINED/partial/plot/"$ALIGNER"/substitutions.0.png $STATS_VCF_FILTERED_COMBINED/partial/plot/"$ALIGNER"_combined_snp_substitutions.png 2>/dev/null
           #mv $STATS_VCF_FILTERED_COMBINED/partial/plot/"$ALIGNER"/tstv_by_qual.0.png $STATS_VCF_FILTERED_COMBINED/partial/plot/"$ALIGNER"_combined_snp_tstv_by_qual.png 2>/dev/null
           #rm -r $STATS_VCF_FILTERED_COMBINED/partial/plot/"$ALIGNER" 2>/dev/null
           #wait
      done
      wait
   fi
   wait
   
   #
   #
   #
      
   if echo $COMBINE_VCF | grep -q "full"; then
   
      for ALIGNER in $ALIGNER_PROGRAM; do
         export ALIGNER=$ALIGNER
         echo "##### $ALIGNER #####"
         source "$IVDP_DIR"/program/01.general_steps/folders.sh
              
           ### MEAN SITE DEPTH
           #vcftools --gzvcf $VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp.vcf.gz \
           #--site-mean-depth --out $STATS_VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp > /dev/null 2>&1
           #wait
           
           
           ### MINOR ALLELE FREQUENCY (MAF)
           #vcftools --gzvcf $VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp.vcf.gz \
           #--freq --out $STATS_VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp > /dev/null 2>&1
           #wait
           
           plink2 --vcf $VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp.vcf.gz \
           --chr-set $CHROM_SET no-xy no-mt \
           --snps-only \
           --min-alleles 2 \
           --max-alleles 2 \
           --set-missing-var-ids @:# \
           --freq \
           --keep-allele-order \
           --const-fid 0 \
           --out $STATS_VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp > /dev/null 2>&1
           wait
           #--freq alt1bins=0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5 \
           rm -r $STATS_VCF_FILTERED_COMBINED/full/*.log 2> /dev/null
           rm -r $STATS_VCF_FILTERED_COMBINED/full/*.nosex 2> /dev/null
           wait
                      
           
           ### MISSING DATA PER INDIVIDUAL
           #vcftools --gzvcf $VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp.vcf.gz \
           #--missing-indv --out $STATS_VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp > /dev/null 2>&1
           #wait
           
           plink2 --vcf $VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp.vcf.gz \
           --chr-set $CHROM_SET no-xy no-mt \
           --snps-only \
           --min-alleles 2 \
           --max-alleles 2 \
           --set-missing-var-ids @:# \
           --missing \
           --keep-allele-order \
           --const-fid 0 \
           --out $STATS_VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp > /dev/null 2>&1
           wait
           rm -r $STATS_VCF_FILTERED_COMBINED/full/*.log 2> /dev/null
           rm -r $STATS_VCF_FILTERED_COMBINED/full/*.nosex 2> /dev/null
           wait    


           ### INBREEDING
            #vcftools --gzvcf $VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp.vcf.gz \
           #--het --out $STATS_VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp > /dev/null 2>&1
           #wait          
           
           plink --vcf $VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp.vcf.gz \
           --chr-set $CHROM_SET no-xy no-mt \
           --snps-only --biallelic-only \
           --set-missing-var-ids @:# \
           --het \
           --keep-allele-order \
           --const-fid 0 \
           --out $STATS_VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp > /dev/null 2>&1
           wait
           rm -r $STATS_VCF_FILTERED_COMBINED/full/*.log 2> /dev/null
           rm -r $STATS_VCF_FILTERED_COMBINED/full/*.nosex 2> /dev/null
           wait
                                 
           for chrom in $(cat $CHROM_NAME); do 
              bcftools stats --threads $THREADS -r $chrom $VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp.vcf.gz > $TMP/fullCombinedVCF_"$ALIGNER"_chrom"$chrom".stats & 
           done
           wait
           rm $TMP/fullCombinedVCF_"$ALIGNER".stats 2> /dev/null
           for chrom in $(cat $CHROM_NAME); do
              echo -e $chrom "\t" snp "\t" $(grep -w "number of SNPs:" $TMP/fullCombinedVCF_"$ALIGNER"_chrom"$chrom".stats | awk -F":" '{print $2}' | sed "s/[[:space:]]\+//g") >> $TMP/fullCombinedVCF_"$ALIGNER".stats
           done
           wait
           unset chrom

           ### STATISTICS
           Rscript -e "inputPath=\"$STATS_VCF_FILTERED_COMBINED/full\";aligner=\"${ALIGNER}\";caller=\"none\";step=\"fullCombinedVCF\";chromNames=\"${CHROM_NAME}\";file1=\"${TMP}/fullCombinedVCF_${ALIGNER}.stats\";file2=\"$STATS_VCF_FILTERED_COMBINED/full/${ALIGNER}_combined_snp.afreq\";file3=\"$STATS_VCF_FILTERED_COMBINED/full/${ALIGNER}_combined_snp.vmiss\";file4=\"$STATS_VCF_FILTERED_COMBINED/full/${ALIGNER}_combined_snp.smiss\";file5=\"$STATS_VCF_FILTERED_COMBINED/full/${ALIGNER}_combined_snp.het\";source(\"$IVDP_DIR/program/01.general_steps/plot_vcf_stats.R\")"
           wait
           bcftools query -f '%CHROM\t%POS\t%CHROM:%POS\n' $VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp.vcf.gz > $VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp.pos
           wait
           Rscript -e "inputFile=\"$VCF_FILTERED_COMBINED/full/${ALIGNER}_combined_snp.pos\";aligner=\"${ALIGNER}\";caller=\"none\";outputPath=\"$STATS_VCF_FILTERED_COMBINED/full\";source(\"$IVDP_DIR/program/01.general_steps/plot_snpDensity.R\")"
           wait
           mv $STATS_VCF_FILTERED_COMBINED/full/Col1.Col0.SNP-Density.plot_DensitySNPs_"$ALIGNER"_none.pdf $STATS_VCF_FILTERED_COMBINED/full/plot_DensitySNPs_"$ALIGNER"_"$CALLER"_combined_snp.pdf
           wait

           gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dAutoRotatePages=/None -sOutputFile=$REPORT_VCF_FILTERED_COMBINED/full/plot_stats_"$ALIGNER"_combined_snp.pdf $STATS_VCF_FILTERED_COMBINED/full/plot_stats_"$ALIGNER"_"$CALLER"_page1.pdf $STATS_VCF_FILTERED_COMBINED/full/plot_DensitySNPs_"$ALIGNER"_"$CALLER"_combined_snp.pdf $STATS_VCF_FILTERED_COMBINED/full/plot_stats_"$ALIGNER"_"$CALLER"_page2.pdf $STATS_VCF_FILTERED_COMBINED/full/plot_stats_"$ALIGNER"_"$CALLER"_page3.pdf
           wait
           rm $STATS_VCF_FILTERED_COMBINED/full/plot_DensitySNPs_"$ALIGNER"_"$CALLER"_combined_snp.pdf $STATS_VCF_FILTERED_COMBINED/full/plot_stats_"$ALIGNER"_"$CALLER"_page* $STATS_VCF_FILTERED_COMBINED/full/Rplots.pdf $VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp.pos 2> /dev/null
           wait
           
           bcftools stats $VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp.vcf.gz > $STATS_VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp.stats
           wait
           
           #plot-vcfstats -P -p $STATS_VCF_FILTERED_COMBINED/full/plot/"$ALIGNER" -t ${ALIGNER}_combined_snp $STATS_VCF_FILTERED_COMBINED/full/"$ALIGNER"_combined_snp.stats 2>/dev/null
           #mv $STATS_VCF_FILTERED_COMBINED/full/plot/"$ALIGNER"/depth.0.png $STATS_VCF_FILTERED_COMBINED/full/plot/"$ALIGNER"_combined_snp_depth.png 2>/dev/null
           #mv $STATS_VCF_FILTERED_COMBINED/full/plot/"$ALIGNER"/substitutions.0.png $STATS_VCF_FILTERED_COMBINED/full/plot/"$ALIGNER"_combined_snp_substitutions.png 2>/dev/null
           #mv $STATS_VCF_FILTERED_COMBINED/full/plot/"$ALIGNER"/tstv_by_qual.0.png $STATS_VCF_FILTERED_COMBINED/full/plot/"$ALIGNER"_combined_snp_tstv_by_qual.png 2>/dev/null
           #rm -r $STATS_VCF_FILTERED_COMBINED/full/plot/"$ALIGNER" 2>/dev/null
           #wait
      done
      wait
   fi
   wait

fi
wait


###############################
##### MULTIQC COMBINE VCF #####
###############################

if [[ "$N_CALLERS" -gt "1" ]]; then
   echo
   echo "MULTIQC REPORT"
   
   if echo $COMBINE_VCF | grep -q "partial"; then
      multiqc -f -n report_vcf_combined_snp.html $STATS_VCF_FILTERED_COMBINED/partial/* -o $REPORT_VCF_FILTERED_COMBINED/partial > /dev/null 2>&1
   fi
   wait
   
   #
   
   if echo $COMBINE_VCF | grep -q "full"; then
      multiqc -f -n report_vcf_combined_snp.html $STATS_VCF_FILTERED_COMBINED/full/* -o $REPORT_VCF_FILTERED_COMBINED/full > /dev/null 2>&1
   fi
   wait

fi
wait


END2=$(date +%s)
DIFF2=$(( $END2 - $START2 ))

echo
echo "#@ COMBINE VCF TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF2/3600)) $(($DIFF2%3600/60)) $(($DIFF2%60)))"
echo "#@#############################################################"