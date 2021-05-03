#!/bin/bash

echo
echo
echo "#@#############################################################"
echo "#@	  DOWNSAMPLING BAM FILES: $OUTPUT_NAME "
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo

START1=$(date +%s)

##### Check Files from the Previous Step
for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
   source "$IVDP_DIR"/program/01.general_steps/folders.sh

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


##### Run DownSampling BAM
n=0
for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
   source "$IVDP_DIR"/program/01.general_steps/folders.sh
   
   mkdir -p "$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/downsampling_bam_"$DOWNSAMPLING_BAM"x/"$ALIGNER" > /dev/null 2>&1
   DOWNSAMPLING_FOLDER="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/downsampling_bam_"$DOWNSAMPLING_BAM"x/"$ALIGNER"
	
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      export i=$( echo $sample | awk -F\\ '{print $1}')
      export SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      export ID_PU=$( echo $sample |awk -F\\ '{print $3}')

      echo "$ALIGNER SAMPLE: $i"
            
      # Depth of Coverage
      #N=$(grep ${i} $STATS_ALIGN/align_stats.txt | awk '{printf "%.0f\n", $2}')
      
      # Genome Size
      G=$(samtools idxstats $ALIGN/${i}_trim_"$ALIGNER"_sorted_RG.bam | cut -f2 | awk 'BEGIN {total=0} {total += $1} END {print total}') && \
      
      # Number of mapped reads
      N_READS=$(samtools idxstats $ALIGN/${i}_trim_"$ALIGNER"_sorted_RG.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {print total}') && \
      
      # Median of read length
      L=$(samtools view -@ $THREADS $ALIGN/${i}_trim_"$ALIGNER"_sorted_RG.bam | awk '{print length($10)}' | datamash median 1) && \
      
      # Number of reads to be sample to achieve a certain coverage (newCoverage * G / L)
      N_READS_SAMPLING=$(echo "scale=4; $DOWNSAMPLING_BAM * $G / $L" | bc) && \
      
      # Proportion of reads to be sampled
      PROPORTION_OF_READS=$(echo "scale=4; $N_READS_SAMPLING / $N_READS" | bc) && \
      
      seed=$(echo $RANDOM) && \
      sambamba view -h -t 8 -s $PROPORTION_OF_READS -f bam --subsampling-seed=$seed -o $DOWNSAMPLING_FOLDER/${i}_trim_"$ALIGNER"_sorted_RG.bam $ALIGN/${i}_trim_"$ALIGNER"_sorted_RG.bam && \
      wait
      
      # Statistics downsampled bam file
      mkdir -p $DOWNSAMPLING_FOLDER/stats > /dev/null 2>&1 
      bedtools genomecov -bga -ibam $DOWNSAMPLING_FOLDER/${i}_trim_"$ALIGNER"_sorted_RG.bam > $DOWNSAMPLING_FOLDER/stats/${i}.coverageBed && \
      wait
      totalBasesGenome=$( samtools view -H $DOWNSAMPLING_FOLDER/${i}_trim_"$ALIGNER"_sorted_RG.bam | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}' ) && \
      wait
      awk '{numberPositions = $3 - $2; basesCovered = numberPositions * $4; print numberPositions, basesCovered}' $DOWNSAMPLING_FOLDER/stats/${i}.coverageBed | awk -v totalBasesGenome=$totalBasesGenome '$2 > 0 {sumCoverage+=$2; totalBasesCovered+=$1} END {printf "%.2f\t%.2f\t%.2f\n", sumCoverage/totalBasesGenome, (totalBasesCovered/totalBasesGenome)*100, sumCoverage/totalBasesCovered}' >> $DOWNSAMPLING_FOLDER/stats/"${i}".cov && \
      wait
      samtools idxstats $DOWNSAMPLING_FOLDER/${i}_trim_"$ALIGNER"_sorted_RG.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {print total}' > $DOWNSAMPLING_FOLDER/stats/"${i}".reads &
      wait
      
      # limit jobs
      if (( $(($((++n)) % $BATCH)) == 0 )) ; then
      wait # wait until all have finished (not optimal, but most times good enough)
      echo $n Files completed
      fi
   done
done
wait


##### Merge statistics Downsampled BAM files

echo "Running statistics of Downsampled BAM files"

for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
   
   DOWNSAMPLING_FOLDER="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/downsampling_bam_"$DOWNSAMPLING_BAM"x/"$ALIGNER"
   touch $DOWNSAMPLING_FOLDER/stats/align_stats_"$DOWNSAMPLING_BAM"x.txt
   echo $'Sample\tDepth\tBreadth_Percentage\tDepth_only_cov_regions\tTotal_Mapped_Reads' > $DOWNSAMPLING_FOLDER/stats/align_stats_"$DOWNSAMPLING_BAM"x.txt
   wait

   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      export i=$( echo $sample | awk -F\\ '{print $1}')
      export SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      export ID_PU=$( echo $sample |awk -F\\ '{print $3}')
   
      paste <(echo "$i") <(echo "$(cat $DOWNSAMPLING_FOLDER/stats/"${i}".cov)") <(echo "$(cat $DOWNSAMPLING_FOLDER/stats/"${i}".reads)") -d '\t' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' >> $DOWNSAMPLING_FOLDER/stats/align_stats_"$DOWNSAMPLING_BAM"x.txt
   
   done
done
wait

rm $DOWNSAMPLING_FOLDER/stats/*.coverageBed > /dev/null 2>&1
rm $DOWNSAMPLING_FOLDER/stats/*.cov > /dev/null 2>&1
rm $DOWNSAMPLING_FOLDER/stats/*.reads > /dev/null 2>&1 
wait


END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ DOWNSAMPLING BAM FILES TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@#############################################################"

#http://crazyhottommy.blogspot.com/2016/05/downsampling-for-bam-files-to-certain.html
#https://github.com/pughlab/wadingpool/wiki/Downsample-bam-file
#https://www.biostars.org/p/5165/
#https://www.biostars.org/p/104063/
#https://bioinformatics.stackexchange.com/questions/13722/calculating-average-coverage-for-bam-files-sequence-data
#https://www.biostars.org/p/356937/
#https://www.ecseq.com/support/ngs/how-to-calculate-the-coverage-for-a-sequencing-experiment
