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
   if [[ "$(ls *sorted_RG.bam 2>/dev/null | wc -l)" -eq "0" ]]; then
      echo "#@ ERROR: No sorted bam files from alignment step in $ALIGN "
      echo "Check if STEP_ALIGNMENT=yes and the log files in $LOG_ALIGN"
      echo "Aborting analysis"
      exit 1
   else
      for sample in $(ls *sorted_RG.bam); do
         #if [[ "$(du $sample | awk '{print $1}')" -lt "1000" ]]; then
         if [[ "$(samtools view $sample | head -n 1000 | wc -l)" == 0 ]]; then
            echo "#@ ERROR: Malformed sorted bam files from alignment step in $ALIGN: $sample"
            echo "Check if STEP_ALIGNMENT=yes and the log files in $LOG_ALIGN"
            echo "Aborting analysis"
            exit 1
         fi
      done
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
n=0
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
   
   
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      export i=$( echo $sample | awk -F\\ '{print $1}')
      export SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      export ID_PU=$( echo $sample |awk -F\\ '{print $3}')

      echo "$ALIGNER SAMPLE: $i"
      source "$IVDP_DIR"/program/01.general_steps/folders.sh


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
      fi   &
      
      
      # limit jobs
      if (( $(($((++n)) % $BATCH)) == 0 )) ; then
      wait # wait until all have finished (not optimal, but most times good enough)
      echo $n Files completed
      fi
  
   done
done
wait
echo


END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ MDUP and SPLITNCIGAR TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@#############################################################"
