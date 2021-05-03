#!/bin/bash

echo
echo
echo "#@#############################################################"
echo "#@	  SUMMARY STATISTICS FASTQ FILES: $OUTPUT_NAME "
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo

START1=$(date +%s)

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
         echo "Check the samples in $INPUT_DIR"
         echo "Aborting analysis"
         exit 1
      fi     
   fi
fi

##### Check Output Folders for This Step
source "$IVDP_DIR"/program/01.general_steps/folders.sh
mkdir -p "$REPORT_TRIM" > /dev/null 2>&1
wait


##### Run QC and Adapter Trimming
echo "##### Running summary statistics of Fastq files"


if [ ! -f "$REPORT_TRIM"/fastq_stats.txt ]; then
   touch "$REPORT_TRIM"/fastq_stats.txt
   echo $'Sample\tTotal Reads\tTotal Bases\tN Bases\tQ20\tQ30\tGC' > "$REPORT_TRIM"/fastq_stats.txt
fi
wait


if [[ "$TYPE" == "se" ]]; then
   n=0
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      export i=$( echo $sample | awk -F\\ '{print $1}')
      export SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      export ID_PU=$( echo $sample |awk -F\\ '{print $3}')
      echo "SAMPLE: $i"
      
      if ! grep -q $i "$REPORT_TRIM"/fastq_stats.txt; then
         paste -d"\t" <(echo $i) <(pigz -fdc ${i}$EXTENSION_SE | $IVDP_DIR/program/03.qc_trim/FastqCount/FastqCount - | tail -n +2) >> "$REPORT_TRIM"/fastq_stats.txt
      fi &
            
      # limit jobs
      if (( $(($((++n)) % 30)) == 0 )) ; then
      wait # wait until all have finished (not optimal, but most times good enough)
      echo $n Files completed
      fi
      
   done
   wait

else
   n=0
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      export i=$( echo $sample | awk -F\\ '{print $1}')
      export SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      export ID_PU=$( echo $sample |awk -F\\ '{print $3}')
      echo "SAMPLE: $i"
      
      if ! grep -q $i "$REPORT_TRIM"/fastq_stats.txt; then
         paste -d"\t" <(echo $i) <(pigz -fdc ${i}$EXTENSION_PE1 | $IVDP_DIR/program/03.qc_trim/FastqCount/FastqCount - | tail -n +2) >> "$REPORT_TRIM"/fastq_stats.txt
      fi &
      
      # limit jobs
      if (( $(($((++n)) % 30)) == 0 )) ; then
      wait # wait until all have finished (not optimal, but most times good enough)
      echo $n Files completed
      fi
      
   done
   wait
fi
wait

Rscript -e "inputPath=\"${REPORT_TRIM}\";source(\"$IVDP_DIR/program/03.qc_trim/fastq_stats_table.R\")"
wait


END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ SUMMARY STATISTICS FASTQ FILES TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@#############################################################"
