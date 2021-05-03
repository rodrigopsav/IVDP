#!/bin/bash

echo
echo
echo "#@#############################################################"
echo "#@ MDUP and SPLITNCIGAR: $OUTPUT_NAME "
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


   if [[ $DOWNSAMPLING_BAM != 0 ]]; then
      export ALIGN="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/downsampling_bam_"$DOWNSAMPLING_BAM"x/"$ALIGNER"
   else
      if [[ -d "$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/downsampling_bam_"$DOWNSAMPLING_BAM"x/"$ALIGNER" ]]; then
         export ALIGN="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/downsampling_bam_"$DOWNSAMPLING_BAM"x/"$ALIGNER"
      else
         export ALIGN=$ALIGN
      fi
   fi
   wait
   
      
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
   source "$IVDP_DIR"/program/01.general_steps/folders.sh

   mkdir -p "$ALIGN_MDUP" > /dev/null 2>&1
   mkdir -p "$LOG_ALIGN_MDUP" > /dev/null 2>&1

done
wait


##### Run Alignment Correction
for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER   
   echo "$ALIGNER SAMPLE: $i"
   source "$IVDP_DIR"/program/01.general_steps/folders.sh


   if [[ $DOWNSAMPLING_BAM != 0 ]]; then
      export ALIGN="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/downsampling_bam_"$DOWNSAMPLING_BAM"x/"$ALIGNER"
   else
      if [[ -d "$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/downsampling_bam_"$DOWNSAMPLING_BAM"x/"$ALIGNER" ]]; then
         export ALIGN="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/downsampling_bam_"$DOWNSAMPLING_BAM"x/"$ALIGNER"
      else
         export ALIGN=$ALIGN
      fi
   fi
   wait
   
   
##### MARK DUPLICATED READS
   if [[ "$MDUP_PROGRAM" == "markdupspark" ]]; then
   	$IVDP_DIR/program/06.mdup_bqsr/1.mdup_markdupspark.sh > $LOG_ALIGN_MDUP/"$i"_mdup.txt 2>&1
   
   elif [[ "$MDUP_PROGRAM" == "picard" ]]; then
   	$IVDP_DIR/program/06.mdup_bqsr/1.mdup_picard.sh > $LOG_ALIGN_MDUP/"$i"_mdup.txt 2>&1
   
   elif [[ "$MDUP_PROGRAM" == "sambamba" ]]; then
   	$IVDP_DIR/program/06.mdup_bqsr/1.mdup_sambamba.sh > $LOG_ALIGN_MDUP/"$i"_mdup.txt 2>&1
   
   else echo ""
   
   fi &&
 
 
##### SPLIT "N" CIGAR

   if [[ "$ANALYSIS" == "1" ||  "$ALIGNER" == "bbmap" || "$ALIGNER" == "bowtie" ]]; then
      echo " "
   else
      $IVDP_DIR/program/06.mdup_bqsr/2.splitncigar.sh > $LOG_ALIGN_MDUP/"$i"_splitncigar.txt 2>&1
   fi &

done
wait

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ MDUP and SPLITNCIGAR TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@#############################################################"
