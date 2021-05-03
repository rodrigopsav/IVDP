#!/bin/bash

echo
echo
echo "#@#############################################################"
echo "#@                 GENE_COUNT: $OUTPUT_NAME "
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo

START1=$(date +%s)

######################
##### GENE COUNT #####
######################

##### Check Files from the Previous Step
for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
   source "$IVDP_DIR"/program/01.general_steps/folders.sh
   
   cd $ALIGN
   if [[ "$(ls ${i}_trim_${ALIGNER}_sorted_RG.bam 2>/dev/null | wc -l)" -eq "0" ]]; then
      echo "#@ ERROR: No sorted bam files from alignment step in $ALIGN: ${i}_trim_${ALIGNER}_sorted_RG.bam"
      echo "Check if STEP_ALIGNMENT=yes and the log files in $LOG_ALIGN"
      echo "Aborting analysis"
      exit 1
   else
      #if [[ "$(du ${i}_trim_${ALIGNER}_sorted_RG.bam | awk '{print $1}')" -lt "1000" ]]; then
      if [[ "$(samtools view ${i}_trim_${ALIGNER}_sorted_RG.bam | head -n 1000 | wc -l)" == 0 ]]; then
         echo "#@ ERROR: Malformed sorted bam files from alignment step in $ALIGN: ${i}_trim_${ALIGNER}_sorted_RG.bam"
         echo "Check if STEP_ALIGNMENT=yes and the log files in $LOG_ALIGN"
         echo "Aborting analysis"
         exit 1
      fi
   fi
done
wait


##### Check Output Folders for This Step
for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
   for GENE_COUNT in $GENECOUNT_PROGRAM; do
      export GENE_COUNT=$GENE_COUNT
      for FEATURE in $FEATURE_TYPE; do
         export FEATURE=$FEATURE
      
         source "$IVDP_DIR"/program/01.general_steps/folders.sh
   
         mkdir -p "$GENECOUNT_DATA" > /dev/null 2>&1
         mkdir -p "$LOG_GENECOUNT" > /dev/null 2>&1
         mkdir -p "$REPORT_GENECOUNT" > /dev/null 2>&1
        
      done
   done      
done
wait


##### Run Gene Count for RNAseq
for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
   for GENE_COUNT in $GENECOUNT_PROGRAM; do
      export GENE_COUNT=$GENE_COUNT
      for FEATURE in $FEATURE_TYPE; do
         export FEATURE=$FEATURE
         echo "$ALIGNER $GENE_COUNT $FEATURE: SAMPLE $i"
         source "$IVDP_DIR"/program/01.general_steps/folders.sh
   
         ##### GENE COUNTS
		 $IVDP_DIR/program/05.gene_count/genecount_"$GENE_COUNT".sh >> $LOG_GENECOUNT/"$i".txt 2>&1 &

      done
   done
done
wait

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ GENE COUNT TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@#############################################################"
