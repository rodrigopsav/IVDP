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
      if [[ "$(ls ${i}_trim_${ALIGNER}_sorted_RG_mdup.bam 2>/dev/null | wc -l)" -eq "0" ]]; then
         echo "#@ ERROR: No mdup.bam files from mark duplicated read step in $ALIGN_MDUP: ${i}_trim_${ALIGNER}_sorted_RG_mdup.bam"
         echo "Check if STEP_MDUP=yes and the log mdup.txt files in $LOG_ALIGN_MDUP"
         echo "Aborting analysis"
         exit 1
      else
         #if [[ "$(du ${i}_trim_${ALIGNER}_sorted_RG_mdup.bam | awk '{print $1}')" -lt "1000" ]]; then
         if [[ "$(samtools view ${i}_trim_${ALIGNER}_sorted_RG_mdup.bam | head -n 1000 | wc -l)" == 0 ]]; then
            echo "#@ ERROR: Malformed mdup.bam files from mark duplicated read step in $ALIGN_MDUP: ${i}_trim_${ALIGNER}_sorted_RG_mdup.bam"
            echo "Check if STEP_MDUP=yes and the log mdup.txt files in $LOG_ALIGN_MDUP"
            echo "Aborting analysis"
            exit 1
         fi
      fi
      
   else
   
      cd $ALIGN_MDUP
      if [[ "$(ls ${i}_trim_${ALIGNER}_sorted_RG_mdup_splitncigar.bam 2>/dev/null | wc -l)" -eq "0" ]]; then
         echo "#@ ERROR: No mdup_splitncigar.bam files from mark duplicated read and splitncigar steps in $ALIGN_MDUP: ${i}_trim_${ALIGNER}_sorted_RG_mdup_splitncigar.bam "
         echo "Check if STEP_MDUP=yes and the log mdup.txt and splitncigar.txt files in $LOG_ALIGN_MDUP"
         echo "Aborting analysis"
         exit 1
      else
         #if [[ "$(du ${i}_trim_${ALIGNER}_sorted_RG_mdup_splitncigar.bam | awk '{print $1}')" -lt "1000" ]]; then
         if [[ "$(samtools view ${i}_trim_${ALIGNER}_sorted_RG_mdup_splitncigar.bam | head -n 1000 | wc -l)" == 0 ]]; then
            echo "#@ ERROR: Malformed mdup.bam files from mark duplicated read and splitncigar steps in $ALIGN_MDUP: ${i}_trim_${ALIGNER}_sorted_RG_mdup_splitncigar.bam"
            echo "Check if STEP_MDUP=yes and the log mdup.txt and splitncigar.txt files in $LOG_ALIGN_MDUP"
            echo "Aborting analysis"
            exit 1
         fi
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
for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
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

done
wait



######################################
##### Check Files from BQSR Step #####
######################################
for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
   source "$IVDP_DIR"/program/01.general_steps/folders.sh

   cd $ALIGN_BQSR
   if [[ "$(ls ${i}_trim_${ALIGNER}_*bqsr.bam 2>/dev/null | wc -l)" -eq "0" ]]; then
      echo "#@ ERROR: No bqsr.bam files from base quality score recalibration step in $ALIGN_BQSR"
      echo "Check if STEP_BQSR=yes and the log bqsr.txt files in $LOG_ALIGN_BQSR"
      echo "Aborting analysis"
      exit 1
   else
      #if [[ "$(du ${i}*bqsr.bam | awk '{print $1}')" -lt "1000" ]]; then
      if [[ "$(samtools view ${i}_trim_${ALIGNER}_*bqsr.bam | head -n 1000 | wc -l)" == 0 ]]; then
         echo "#@ ERROR: Malformed bqsr.bam files from base quality score recalibration step in $ALIGN_BQSR: ${i}_trim_${ALIGNER}_*bqsr.bam"
         echo "Check if STEP_BQSR=yes and the log bqsr.txt files in $LOG_ALIGN_BQSR"
         echo "Aborting analysis"
         exit 1
      fi
   fi
done
wait



##################################################
echo "#@##### FILES FOR VARIANT CALLING"
##################################################

for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
   echo "#@######################################"
   echo "#@ $ALIGNER"
   source $IVDP_DIR/program/01.general_steps/folders.sh

   splitncigar="$(ls $ALIGN_MDUP/*splitncigar.bam 2> /dev/null | wc -l)"
   mdup="$(ls $ALIGN_MDUP/*RG_mdup.bam 2> /dev/null | wc -l)"
   bqsr="$(ls $ALIGN_BQSR/*bqsr.bam 2> /dev/null | wc -l)"

##### MDUP AND BQSR FILES
   if [[ "$splitncigar" != 0 ]]; then
      rm -r $ALIGN_MDUP/${i}_trim_"$ALIGNER"_sorted_RG_mdup.bam 2> /dev/null
      rm -r $ALIGN_MDUP/${i}_trim_"$ALIGNER"_sorted_RG_mdup.bam.bai 2> /dev/null
      mv $ALIGN_MDUP/${i}_trim_"$ALIGNER"_sorted_RG_mdup_splitncigar.bam $ALIGN_MDUP/${i}_trim_"$ALIGNER"_sorted_RG_mdup.bam
      mv $ALIGN_MDUP/${i}_trim_"$ALIGNER"_sorted_RG_mdup_splitncigar.bam.bai $ALIGN_MDUP/${i}_trim_"$ALIGNER"_sorted_RG_mdup.bam.bai
   fi
   wait


   if [[ "$mdup" != 0 ]] && [[ "$bqsr" != 0 ]]; then
      echo "#@ MDUP, AND BQSR files exist"
      echo "#@ Use $ALIGN_BQSR/*bqsr.bam for variant calling (except for FREEBAYES)"
      echo "#@ Use $ALIGN_MDUP/*mdup.bam for variant calling with FREEBAYES"
        
        
   ##### RG FILES
   elif [[ "$mdup" != 0 ]] && [[ "$bqsr" == 0 ]]; then
      echo "#@ ONLY MDUP files exist"
      echo "#@ Use $ALIGN_MDUP/*mdup.bam for variant calling"
   
   else
      echo "#@ ERROR: THERE ARE NO MDUP OR BQSR BAM FILES"
      echo "#@ PLEASE CHECK IF THE BAM FILES FROM ALIGNMENT STEP ARE OK (*RG.bam files in $ALIGN SUBFOLDERS)"
      echo "#@ PLEASE CHECK THE LOG FILES IN:"
      echo "$LOG_ALIGN"
      echo "$LOG_ALIGN_MDUP"
      echo "$LOG_ALIGN_BQSR"
      exit 1
   fi

echo "#@######################################"
echo
done
wait

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ BQSR TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@#############################################################"
