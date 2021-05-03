#!/bin/bash

echo
echo
echo "#@#############################################################"
echo "#@ BQSR: $OUTPUT_NAME "
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo

START1=$(date +%s)

################################################
##### ALIGNMENT CORERECTION OF FASTQ FILES #####
################################################


##### Check Files from the Previous Step
for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
   source "$IVDP_DIR"/program/01.general_steps/folders.sh

   if [[ "$ANALYSIS" == "1" ]]; then
      cd $ALIGN_MDUP
      if [[ "$(ls *mdup.bam 2>/dev/null | wc -l)" -eq "0" ]]; then
         echo "#@ ERROR: No mdup.bam files from mark duplicated read step in $ALIGN_MDUP "
         echo "Check if STEP_MDUP=yes and the log mdup.txt files in $LOG_ALIGN_MDUP"
         echo "Aborting analysis"
         exit 1
      else
         for sample in $(ls *mdup.bam); do
            #if [[ "$(du $sample | awk '{print $1}')" -lt "1000" ]]; then
            if [[ "$(samtools view $sample | head -n 1000 | wc -l)" == 0 ]]; then 
               echo "#@ ERROR: Malformed mdup.bam files from mark duplicated read step in $ALIGN_MDUP: $sample"
               echo "Check if STEP_MDUP=yes and the log mdup.txt files in $LOG_ALIGN_MDUP"
               echo "Aborting analysis"
               exit 1
            fi
         done
      fi
      
   else
   
      cd $ALIGN_MDUP
      if [[ "$(ls *splitncigar.bam 2>/dev/null | wc -l)" -eq "0" ]]; then
         echo "#@ ERROR: No mdup_splitncigar.bam files from mark duplicated read and splitncigar steps in $ALIGN_MDUP "
         echo "Check if STEP_MDUP=yes and the log mdup.txt and splitncigar.txt files in $LOG_ALIGN_MDUP"
         echo "Aborting analysis"
         exit 1
      else
         for sample in $(ls *splitncigar.bam); do
            #if [[ "$(du $sample | awk '{print $1}')" -lt "1000" ]]; then
            if [[ "$(samtools view $sample | head -n 1000 | wc -l)" == 0 ]]; then
               echo "#@ ERROR: Malformed mdup.bam files from mark duplicated read and splitncigar steps in $ALIGN_MDUP: $sample"
               echo "Check if STEP_MDUP=yes and the log mdup.txt and splitncigar.txt files in $LOG_ALIGN_MDUP"
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

   mkdir -p "$ALIGN_BQSR" > /dev/null 2>&1
   mkdir -p "$LOG_ALIGN_BQSR" > /dev/null 2>&1

done
wait


##### Run Alignment Correction
n=0
for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
	
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      export i=$( echo $sample | awk -F\\ '{print $1}')
      export SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      export ID_PU=$( echo $sample |awk -F\\ '{print $3}')
      
      echo "$ALIGNER SAMPLE: $i"
      source "$IVDP_DIR"/program/01.general_steps/folders.sh


##### BASE QUALITY SCORE RECALIBRATION

      if [[ "$VARIANTS" == "none" || "$ALIGNER" == "bbmap" || "$ALIGNER" == "bowtie" ]]; then
         echo " "
      else
         if [[ "$ANALYSIS" == "1" && "$BQSR_PROGRAM" == "bqsr" ]]; then
              $IVDP_DIR/program/06.mdup_bqsr/3.bqsr_wgs.sh > $LOG_ALIGN_BQSR/"$i"_bqsr.txt 2>&1

         elif [[ "$ANALYSIS" == "1" && "$BQSR_PROGRAM" == "bqsrgatk3" ]]; then
              $IVDP_DIR/program/06.mdup_bqsr/3.bqsr_gatk3_wgs.sh > $LOG_ALIGN_BQSR/"$i"_bqsr.txt 2>&1
                       
         elif [[ "$ANALYSIS" == "1" && "$BQSR_PROGRAM" == "bqsrspark" ]]; then
              $IVDP_DIR/program/06.mdup_bqsr/3.bqsr_wgs_spark.sh > $LOG_ALIGN_BQSR/"$i"_bqsr.txt 2>&1
         
         elif [[ "$ANALYSIS" == "2" && "$BQSR_PROGRAM" == "bqsr" ]]; then
              $IVDP_DIR/program/06.mdup_bqsr/3.bqsr_rnaseq.sh > $LOG_ALIGN_BQSR/"$i"_bqsr.txt 2>&1
         
         elif [[ "$ANALYSIS" == "2" && "$BQSR_PROGRAM" == "bqsrspark" ]]; then
              $IVDP_DIR/program/06.mdup_bqsr/3.bqsr_rnaseq_spark.sh > $LOG_ALIGN_BQSR/"$i"_bqsr.txt 2>&1
         
         else echo " "
         fi
      fi &


      # limit jobs
      if (( $(($((++n)) % $BATCH)) == 0 )) ; then
      wait # wait until all have finished (not optimal, but most times good enough)
      echo $n Files completed
      
      fi

   done
done
wait
echo


##################################################
echo "#@##### FILES FOR VARIANT CALLING"
##################################################
###################
for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
   echo "#@ $ALIGNER"
   source $IVDP_DIR/program/01.general_steps/folders.sh

   $IVDP_DIR/program/06.mdup_bqsr/4.final_bamfiles.sh

   echo "#@######################################"
   echo
   
done
wait

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ BQSR TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@#############################################################"
