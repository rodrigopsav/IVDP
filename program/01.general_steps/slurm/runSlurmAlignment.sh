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
   if [[ "$(ls ${i}*fastq* 2>/dev/null | wc -l)" -eq "0" ]]; then
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
   
      if [[ "$(ls ${i}_trim.fastq.gz 2>/dev/null | wc -l)" -eq "0" ]]; then
         echo "#@ ERROR: No fastq files from adapter trimming step in $TRIM: ${i}_trim.fastq.gz"
         echo "Check if STEP_QC_TRIM=yes and if INPUT_DIR is a valid path in parameter file"
         echo "Check the log files in $LOG_TRIM"
         echo "Aborting analysis"
         exit 1
      else
         #if [[ "$(du ${i}_trim.fastq.gz | awk '{print $1}')" -lt "1000" ]]; then
         if [[ "$(zcat ${i}_trim.fastq.gz | head -n 1000 | wc -l)" == 0 ]]; then
            echo "#@ ERROR: Malformed fastq files from adapter trimming step in $TRIM: ${i}_trim.fastq.gz"
            echo "Check if STEP_QC_TRIM=yes and if INPUT_DIR is a valid path in parameter file"
            echo "Check the log files in $LOG_TRIM"
            echo "Aborting analysis"
            exit 1
         fi
      fi
      
   else
   
      if [[ "$(ls ${i}.*_paired.fastq.gz 2>/dev/null | wc -l)" -eq "0" ]]; then
         echo "#@ ERROR: No fastq files from adapter trimming step in $TRIM: ${i}"
         echo "Check if STEP_QC_TRIM=yes and if INPUT_DIR is a valid path in parameter file"
         echo "Check the log files in $LOG_TRIM"
         echo "Aborting analysis"
         exit 1
      else
         #if [[ "$(du ${i}.R1_paired.fastq.gz | awk '{print $1}')" -lt "1000" && "$(du ${i}.R2_paired.fastq.gz | awk '{print $1}')" -lt "1000" ]]; then
         if [[ "$(zcat ${i}.R1_paired.fastq.gz | head -n 1000 | wc -l)" == 0 && "$(zcat ${i}.R2_paired.fastq.gz | head -n 1000 | wc -l)" == 0 ]]; then
            echo "#@ ERROR: Malformed fastq files from adapter trimming step in $TRIM: ${i}.R1_paired.fastq.gz and ${i}.R2_paired.fastq.gz"
            echo "Check if STEP_QC_TRIM=yes and if INPUT_DIR is a valid path in parameter file"
            echo "Check the log files in $LOG_TRIM"
            echo "Aborting analysis"
            exit 1
         fi
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
for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
   export eval "$ALIGNER"_index="$IVDP_DIR"/dataBases/genome_idx/"$ASSEMBLY"_index_"$ALIGNER"
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
      awk '{numberPositions = $3 - $2; basesCovered = numberPositions * $4; print numberPositions, basesCovered}'  $STATS_ALIGN/${i}.coverageBed | awk -v totalBasesGenome=$totalBasesGenome '$2 > 0 {sumCoverage+=$2; totalBasesCovered+=$1} END {printf "%.2f\t%.2f\t%.2f\n", sumCoverage/totalBasesGenome, (totalBasesCovered/totalBasesGenome)*100, sumCoverage/totalBasesCovered}' >> $STATS_ALIGN/"${i}".cov
      wait
      
   fi &
   
done
wait

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ ALIGNMENT TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@#############################################################"
