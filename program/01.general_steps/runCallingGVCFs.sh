#!/bin/bash

echo
echo
echo "#@#############################################################"
echo "#@           VARIANT CALLING GVCF MODE: $OUTPUT_NAME "
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo

START1=$(date +%s)

########################################
##### VARIANT CALLING FILES: GVCFs #####
########################################

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
   
   mkdir -p "$GVCF_RAW"/"$ALIGNER"_"$CALLER" > /dev/null 2>&1
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
      
      for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
         export i=$( echo $sample | awk -F\\ '{print $1}')
         echo "##### $ALIGNER - $CALLER: sample $i #####"
         
         if [[ "$CALLER" == "gatk4CombineGVCFs" || "$CALLER" == "gatk4GenomicsDBImport" ]]; then
            $IVDP_DIR/program/07.variant_caller/gatk4GVCF.sh > $LOG_VCF_RAW/"$ALIGNER"_"$CALLER".txt 2>&1
         else
            $IVDP_DIR/program/07.variant_caller/gatk3GVCF.sh > $LOG_VCF_RAW/"$ALIGNER"_"$CALLER".txt 2>&1         
         fi &
         
         # limit jobs
         if (( $(($((++n)) % $((BATCH + 30)) )) == 0 )) ; then
         wait # wait until all have finished (not optimal, but most times good enough)
         echo $n Files completed
         fi

      done
   done
   wait
         
### CALL_BY_CHROM == yes
else
   
   for ALIGNER in $ALIGNER_PROGRAM; do
      export ALIGNER=$ALIGNER
      source "$IVDP_DIR"/program/01.general_steps/folders.sh 
      export VCF_TMP=$(echo $VCF_RAW)/"$ALIGNER"_"$CALLER"_"${ANALYSIS_ID}"
      mkdir -p $VCF_TMP > /dev/null 2>&1
      n=0

      for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
         export i=$( echo $sample | awk -F\\ '{print $1}')
         echo "##### $ALIGNER - $CALLER: sample $i #####"
         
         for chr in $(cat $CHROM_NAME); do
            export chr=$chr
            
            if [[ "$CALLER" == "gatk4CombineGVCFs" || "$CALLER" == "gatk4GenomicsDBImport" ]]; then
               $IVDP_DIR/program/07.variant_caller/by_chrom/gatk4GVCF_bychrom.sh > $LOG_VCF_RAW/"$ALIGNER"_"$CALLER"_"$i"_"$chr".txt 2>&1
            else
               $IVDP_DIR/program/07.variant_caller/by_chrom/gatk3GVCF_bychrom.sh > $LOG_VCF_RAW/"$ALIGNER"_"$CALLER"_"$i"_"$chr".txt 2>&1
            fi &
            
            # limit jobs
            if (( $(($((++n)) % $((BATCH + 30)) )) == 0 )) ; then
            wait # wait until all have finished (not optimal, but most times good enough)
            echo $n Files completed
            fi
            
         done
      done
   done
         
fi
wait


##### CONCATENATE CHROMOSOMES OF SAME SAMPLE
if [[ "${CALL_BY_CHROM,,}" == 'yes' ]]; then
   for ALIGNER in $ALIGNER_PROGRAM; do
      export ALIGNER=$ALIGNER
      source "$IVDP_DIR"/program/01.general_steps/folders.sh 
      export VCF_TMP=$(echo $VCF_RAW)/"$ALIGNER"_"$CALLER"_"${ANALYSIS_ID}"
      
      for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
         export i=$( echo $sample | awk -F\\ '{print $1}')
   
         # Get head of gvcfs
         grep ^"#" $(ls "${VCF_TMP}"/"${i}"_*.g.vcf | head -n 1) > $GVCF_RAW/"$ALIGNER"_"$CALLER"/"$i".g.vcf
         wait
         
         # Concatenate gvcfs by sample
         for chr in $(cat $CHROM_NAME); do
            export chr=$chr
            grep -v '^#' "${VCF_TMP}"/"${i}"_"${chr}".g.vcf >> $GVCF_RAW/"$ALIGNER"_"$CALLER"/"$i".g.vcf
         done
         wait
      
         # Index GVCF
         gatk IndexFeatureFile -I $GVCF_RAW/"$ALIGNER"_"$CALLER"/"$i".g.vcf > /dev/null 2>&1
         wait

      done
      wait
         
      rm -r $VCF_TMP
      wait
      
   done
fi
wait

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "@ VARIANT CALLING GVCF MODE TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "@#############################################################"
