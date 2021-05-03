#!/bin/bash

source "$IVDP_DIR"/program/01.general_steps/folders.sh
if [[ $STEP_QC_TRIM == "yes" ]]; then
   export RAW_DATA=$TRIM
else
   
   if [[ -d $TRIM ]]; then
      export RAW_DATA=$TRIM
   else
   
      if [[ $DOWNSAMPLING_FASTQ != 0 ]]; then
         export RAW_DATA="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/downsampling_fastq_"$DOWNSAMPLING_FASTQ"x
      else
      
         if [[ -d "$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/downsampling_fastq_"$DOWNSAMPLING_FASTQ"x ]]; then
            export RAW_DATA="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/downsampling_fastq_"$DOWNSAMPLING_FASTQ"x
         else
            export RAW_DATA=$INPUT_DIR
         fi
      fi
   fi
fi
wait


echo
echo
echo "#@#############################################################"
echo "#@                    ALIGNMENT: $OUTPUT_NAME "
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo

START1=$(date +%s)

####################################
##### ALIGNMENT OF FASTQ FILES #####
####################################


### Check if raw data exist
if [[ -d "$RAW_DATA" ]]; then
   cd $RAW_DATA
   if [[ "$(ls *fastq* 2>/dev/null | wc -l)" -eq "0" ]]; then
      echo "#@ ERROR: No fastq files in $RAW_DATA"
      echo "Check if INPUT_DIR is a valid path in parameter file or if $RAW_DATA exits"
      echo "Aborting analysis"
      exit 1
   fi
fi
wait
   

### If trimmed data exist, check:
if [[ -d $TRIM ]]; then
   cd $TRIM
   if [[ "$TYPE" == "se" ]]; then
   
      if [[ "$(ls *trim.fastq.gz 2>/dev/null | wc -l)" -eq "0" ]]; then
         echo "#@ ERROR: No fastq files from adapter trimming step in $TRIM "
         echo "Check if STEP_QC_TRIM=yes and if INPUT_DIR is a valid path in parameter file"
         echo "Check the log files in $LOG_TRIM"
         echo "Aborting analysis"
         exit 1
      else
         for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
         export i=$( echo $sample | awk -F\\ '{print $1}')
         export SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
         export ID_PU=$( echo $sample |awk -F\\ '{print $3}')
         
            #if [[ "$(du ${i}_trim.fastq.gz | awk '{print $1}')" -lt "1000" ]]; then
            if [[ "$(zcat ${i}_trim.fastq.gz | head -n 1000 | wc -l)" == 0 ]]; then
               echo "#@ ERROR: Malformed fastq files from adapter trimming step in $TRIM: ${i}_trim.fastq.gz"
               echo "Check if STEP_QC_TRIM=yes and if INPUT_DIR is a valid path in parameter file"
               echo "Check the log files in $LOG_TRIM"
               echo "Aborting analysis"
               exit 1
            fi
         done
      fi
      
   else
      
      if [[ "$(ls *paired.fastq.gz 2>/dev/null | wc -l)" -eq "0" ]]; then
         echo "#@ ERROR: No fastq files from adapter trimming step in $TRIM "
         echo "Check if STEP_QC_TRIM=yes and if INPUT_DIR is a valid path in parameter file"
         echo "Check the log files in $LOG_TRIM"
         echo "Aborting analysis"
         exit 1
      else
         for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
         export i=$( echo $sample | awk -F\\ '{print $1}')
         export SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
         export ID_PU=$( echo $sample |awk -F\\ '{print $3}')
         
            #if [[ "$(du ${i}.R1_paired.fastq.gz | awk '{print $1}')" -lt "1000" && "$(du ${i}.R2_paired.fastq.gz | awk '{print $1}')" -lt "1000" ]]; then
            if [[ "$(zcat ${i}.R1_paired.fastq.gz | head -n 1000 | wc -l)" == 0 && "$(zcat ${i}.R2_paired.fastq.gz | head -n 1000 | wc -l)" == 0 ]]; then
               echo "#@ ERROR: Malformed fastq files from adapter trimming step in $TRIM: ${i}.R1_paired.fastq.gz and ${i}.R2_paired.fastq.gz"
               echo "Check if STEP_QC_TRIM=yes and if INPUT_DIR is a valid path in parameter file"
               echo "Check the log files in $LOG_TRIM"
               echo "Aborting analysis"
               exit 1
            fi
         done
      fi
   fi
fi
wait


##### Check Output Folders for This Step
for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
   source "$IVDP_DIR"/program/01.general_steps/folders.sh
   
   mkdir -p "$ALIGN" > /dev/null 2>&1
   mkdir -p "$LOG_ALIGN" > /dev/null 2>&1
   mkdir -p "$STATS_ALIGN" > /dev/null 2>&1
   mkdir -p "$REPORT_ALIGN" > /dev/null 2>&1
      
done
wait
 

##### Run Alignment
n=0
for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
	
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      export i=$( echo $sample | awk -F\\ '{print $1}')
      export SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      export ID_PU=$( echo $sample |awk -F\\ '{print $3}')
      echo "$ALIGNER SAMPLE: $i"
      source "$IVDP_DIR"/program/01.general_steps/folders.sh
      

##### ALIGNMENT
      if [[ "$ALIGNER" == "star" ]]; then
         if [[ "${ANNOTATION^^}" == "none" ]]; then
            $IVDP_DIR/program/04.aligner/star_"$TYPE"_without_annotation.sh > $LOG_ALIGN/"$i".txt 2>&1
         else
            $IVDP_DIR/program/04.aligner/star_"$TYPE".sh > $LOG_ALIGN/"$i".txt 2>&1
         fi
      else
		 $IVDP_DIR/program/04.aligner/"$ALIGNER"_"$TYPE".sh > $STATS_ALIGN/"$i".txt 2>&1
         cp $STATS_ALIGN/"$i".txt $LOG_ALIGN/"$i".txt
      fi &&

##### ADD OR REPLACE READ GROUP
      $IVDP_DIR/program/04.aligner/add_read_group.sh >> $LOG_ALIGN/"$i".txt 2>&1 &&


##### ALIGNMENT STATISTICS
      if [[ "$STATS_ALIGNMENT" == "yes" ]]; then
         # Coverage with bedtools (chrom start end cov)
         bedtools genomecov -bga -ibam $ALIGN/${i}_trim_"$ALIGNER"_sorted_RG.bam > $STATS_ALIGN/${i}.coverageBed
         wait
         
         # Coverage with bedtools at each genome position (chrom pos cov); similar samtools depth; final file is much heavier than previous one.
         #bedtools genomecov -d -ibam $ALIGN/${i}_trim_"$ALIGNER"_sorted_RG.bam > $STATS_ALIGN/${i}.coverage 
         #wait
         
         # Proportion of coverage by depth
         #bedtools genomecov -ibam $ALIGN/${i}_trim_"$ALIGNER"_sorted_RG.bam -g $CHROM_SIZE > $STATS_ALIGN/${i}.propcoverage
         bedtools genomecov -ibam $ALIGN/${i}_trim_"$ALIGNER"_sorted_RG.bam > $STATS_ALIGN/${i}.propcoverage
         wait
         
         # Bam Stats with samtools
         samtools stats -@ $THREADS $ALIGN/"${i}"_trim_"$ALIGNER"_sorted_RG.bam > $STATS_ALIGN/"${i}".stats
         wait
         
         # Bam depth and breadth of coverage with samtools (depth of coverage | breadth | depth given only covered positions)
         #totalBasesGenome=$( samtools view -H $ALIGN/"${i}"_trim_"$ALIGNER"_sorted_RG.bam | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}' )
         #samtools depth $ALIGN/"${i}"_trim_"$ALIGNER"_sorted_RG.bam | awk -v totalBasesGenome=$totalBasesGenome '{sumCoverage+=$3; if($3>0) totalBasesCovered+=1} END {printf "%.2f\t%.2f\t%.2f\n", sumCoverage/totalBasesGenome, (totalBasesCovered/totalBasesGenome)*100, sumCoverage/totalBasesCovered}' >> $STATS_ALIGN/"${i}".cov
         #wait
         
         # Bam depth and breadth of coverage with bedtools (depth of coverage | breadth | depth given only covered positions)
         totalBasesGenome=$( samtools view -H $ALIGN/"${i}"_trim_"$ALIGNER"_sorted_RG.bam | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}' )
         awk '{numberPositions = $3 - $2; basesCovered = numberPositions * $4; print numberPositions, basesCovered}' $STATS_ALIGN/${i}.coverageBed | awk -v totalBasesGenome=$totalBasesGenome '$2 > 0 {sumCoverage+=$2; totalBasesCovered+=$1} END {printf "%.2f\t%.2f\t%.2f\n", sumCoverage/totalBasesGenome, (totalBasesCovered/totalBasesGenome)*100, sumCoverage/totalBasesCovered}' >> $STATS_ALIGN/"${i}".cov
         wait
         
      fi &

      # limit jobs
      if (( $(($((++n)) % $BATCH)) == 0 )) ; then
      wait # wait until all have finished (not optimal, but most times good enough)
      echo $n Files completed
      fi

   done
done
wait


#####################################
##### MERGE ALIGMENT STATISTICS #####
#####################################
if [[ "$STATS_ALIGNMENT" == "yes" ]]; then
   for ALIGNER in $ALIGNER_PROGRAM; do
      export ALIGNER=$ALIGNER
	
      for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
         export i=$( echo $sample | awk -F\\ '{print $1}')
         export SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
         export ID_PU=$( echo $sample |awk -F\\ '{print $3}')
         source "$IVDP_DIR"/program/01.general_steps/folders.sh

         $IVDP_DIR/program/04.aligner/align_stats.sh >> $LOG_ALIGN/"$i".txt 2>&1
         wait
         
         Rscript -e "inputPath=\"${REPORT_ALIGN}\";aligner=\"$ALIGNER\";source(\"$IVDP_DIR/program/04.aligner/align_stats_plots_table.R\")"
         wait         
         
      done
   done
fi


#############################
##### MULTIQC ALIGNMENT #####
#############################
echo
echo "MULTIQC REPORT"

for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
   source "$IVDP_DIR"/program/01.general_steps/folders.sh
   
   $IVDP_DIR/program/04.aligner/multiqc_aligner.sh > /dev/null 2>&1 &
   
done
wait

# Delete aligner stats
for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
	
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      export i=$( echo $sample | awk -F\\ '{print $1}')
      export SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      export ID_PU=$( echo $sample |awk -F\\ '{print $3}')
      source "$IVDP_DIR"/program/01.general_steps/folders.sh
      
      rm -r $STATS_ALIGN/"$i".txt 2> /dev/null

   done
done
wait


END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ ALIGNMENT TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@#############################################################"
