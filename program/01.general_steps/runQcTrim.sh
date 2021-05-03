#!/bin/bash

if [[ $DOWNSAMPLING_FASTQ != 0 ]]; then
   export INPUT_DIR="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/downsampling_fastq_"$DOWNSAMPLING_FASTQ"x
else
   if [[ -d "$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/downsampling_fastq_"$DOWNSAMPLING_FASTQ"x ]]; then
      export INPUT_DIR="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/downsampling_fastq_"$DOWNSAMPLING_FASTQ"x
   else
      export INPUT_DIR=$INPUT_DIR
   fi
fi
wait


echo
echo
echo "#@#############################################################"
echo "#@	  QUALITY CONTROL AND TRIMMING: $OUTPUT_NAME "
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo

START1=$(date +%s)

#############################################
####### QUALITY CONTROL AND TRIMMING ########
#############################################


##### Check Files from the Previous Step
source "$IVDP_DIR"/program/01.general_steps/folders.sh
cd $INPUT_DIR
if [[ "$TYPE" == "se" ]]; then
   if [[ "$(ls *$EXTENSION_SE | wc -l)" == 0 ]]; then
      echo "#@ ERROR: No input fastq files in $INPUT_DIR "
      echo "Check the parameter file (INPUT_DIR variable)"
      echo "Aborting analysis"
      exit 1
   fi
else
   if [[ "$(ls *$EXTENSION_PE1 | wc -l)" == 0 ]] && [[ "$(ls *$EXTENSION_PE2 | wc -l)" == 0 ]]; then
      echo "#@ ERROR: No input fastq files in $INPUT_DIR "
      echo "Check the parameter file (INPUT_DIR variable)"
      echo "Aborting analysis"
      exit 1
   else
      if [[ "$(ls *$EXTENSION_PE1 | wc -l)" != "$(ls *$EXTENSION_PE2 | wc -l)" ]]; then
         echo "#@ ERROR: Number of Paired-end1 files does not match with Paired-end files "
         echo "Check the samples in $INPUT_DIR and the parameter file (variables EXTENSION_PE1 and EXTENSION_PE2)"
         echo "Aborting analysis"
         exit 1
      fi     
   fi
fi


##### Check Output Folders for This Step
source "$IVDP_DIR"/program/01.general_steps/folders.sh
mkdir -p "$TMP" 2> /dev/null
mkdir -p "$TRIM" > /dev/null 2>&1
mkdir -p "$LOG_TRIM" > /dev/null 2>&1
mkdir -p "$STATS_QC1" > /dev/null 2>&1
mkdir -p "$STATS_QC2" > /dev/null 2>&1
mkdir -p "$STATS_TRIM" > /dev/null 2>&1
mkdir -p "$REPORT_TRIM" > /dev/null 2>&1
wait


##### Run QC and Adapter Trimming
echo "##### Running QC and Adapter Trimming"
n=0
for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
   export i=$( echo $sample | awk -F\\ '{print $1}')
   export SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
   export ID_PU=$( echo $sample |awk -F\\ '{print $3}')
   echo "SAMPLE: $i"
   
   $IVDP_DIR/program/03.qc_trim/fastqc_raw_"$TYPE".sh > $LOG_TRIM/"$i".txt 2>&1 && \
   $IVDP_DIR/program/03.qc_trim/trim_"$TYPE".sh > $STATS_TRIM/"$i".txt 2>&1 && \
   cat $STATS_TRIM/"$i".txt >> $LOG_TRIM/"$i".txt && \
   $IVDP_DIR/program/03.qc_trim/fastqc_trim_"$TYPE".sh >> $LOG_TRIM/"$i".txt 2>&1 &
   
   # limit jobs
   if (( $(($((++n)) % $BATCH)) == 0 )) ; then
   wait # wait until all have finished (not optimal, but most times good enough)
   echo $n Files completed
   fi

done
wait
echo

##################################################
###### MULTIQC ORIGINAL AND TRIMMED FILES ########
##################################################
echo "MULTIQC REPORT"
$IVDP_DIR/program/03.qc_trim/multiqc_fastq_raw.sh > /dev/null 2>&1 &
$IVDP_DIR/program/03.qc_trim/multiqc_fastq_trimmed.sh > /dev/null 2>&1 &	
$IVDP_DIR/program/03.qc_trim/multiqc_trimmomatic.sh > /dev/null 2>&1 &
wait


END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ QUALITY CONTROL AND TRIMMING TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@#############################################################"
