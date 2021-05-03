#!/bin/bash

echo
echo
echo "#@#############################################################"
echo "#@	  DOWNSAMPLING FASTQ FILES: $OUTPUT_NAME "
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
mkdir -p "$TMP" 2> /dev/null
mkdir -p "$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/downsampling_fastq_"$DOWNSAMPLING_FASTQ"x > /dev/null 2>&1
DOWNSAMPLING_FOLDER="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/downsampling_fastq_"$DOWNSAMPLING_FASTQ"x
wait

# https://bioinformatics.stackexchange.com/questions/935/fast-way-to-count-number-of-reads-and-number-of-bases-in-a-fastq-file
# https://www.biostars.org/p/243552/

# C=(L * N) / G
# C = coverage (DOWNSAMPLING_FASTQ)
# L = read length
# N = read count  (N_READS)
# G = total genome size

# If the coverage is specified, the new read counts will be:
# READ_COUNT_DOWNSAMPLING=(C * G / L)


if [[ "$TYPE" == "se" ]]; then
   n=0
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      export i=$( echo $sample | awk -F\\ '{print $1}') && \
      export SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}') && \
      export ID_PU=$( echo $sample |awk -F\\ '{print $3}') && \
      echo "SAMPLE: $i"
      
      #pigz -fdc ${i}$EXTENSION_SE | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > $TMP/${i}.readCount && \
      pigz -fdc ${i}$EXTENSION_SE | awk '{if(NR%4==2) print length($1)}' > $TMP/${i}.readCount && \
      
      # Genome Size
      G=$(awk '{s+=$2} END {print s}' $REFERENCE.fai) && \
      
      # Number of original reads
      #N_READS=$(cat $TMP/${i}.readCount | awk '{s+=$1} END {print s}') && \
      N_READS=$(cat $TMP/${i}.readCount | wc -l) && \
      
      # Average read length
      #L=$(cat $TMP/${i}.readCount | awk '{s+=$1;m+=$1*$2} END {printf "%.0f\n", m/s}') && \
      L=$(cat $TMP/${i}.readCount | datamash median 1) && \
      
      # Number of sampling reads
      READ_COUNT_DOWNSAMPLING=$(echo "scale=0; $DOWNSAMPLING_FASTQ * $G / $L" |bc) && \
      seed=$(echo $RANDOM) && \
      
      if [[ "$N_READS" -gt "$READ_COUNT_DOWNSAMPLING" ]]; then
         if echo $EXTENSION_SE | grep -q "gz"; then
            NEW_SE=$(echo $EXTENSION_SE | rev | cut -c 4- | rev)
            wait
            seqtk sample -s$seed $INPUT_DIR/"${i}"$EXTENSION_SE $READ_COUNT_DOWNSAMPLING > $DOWNSAMPLING_FOLDER/"$i"$NEW_SE
            wait
            gzip -f $DOWNSAMPLING_FOLDER/"$i"$NEW_SE
            wait
         else
            seqtk sample -s$seed $INPUT_DIR/"${i}"$EXTENSION_SE $READ_COUNT_DOWNSAMPLING > $DOWNSAMPLING_FOLDER/"$i"$EXTENSION_SE
            wait
         fi
      fi && \
      wait
      
      unset READ_COUNT && \
      unset READ_LENGTH && \
      unset GENOME_SIZE && \
      unset READ_COUNT_DOWNSAMPLING && \
      unset seed && \
      rm -r $TMP/${i}.readCount 2> /dev/null &
      wait
      
      # limit jobs
      if (( $(($((++n)) % $BATCH)) == 0 )) ; then
      wait # wait until all have finished (not optimal, but most times good enough)
      echo $n Files completed
      fi
      
   done
   wait

else

   n=0
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      export i=$( echo $sample | awk -F\\ '{print $1}') && \
      export SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}') && \
      export ID_PU=$( echo $sample |awk -F\\ '{print $3}') && \
      echo "SAMPLE: $i" && \

      #pigz -fdc ${i}$EXTENSION_PE1 | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > $TMP/${i}.readCount && \
      pigz -fdc ${i}$EXTENSION_PE1 | awk '{if(NR%4==2) print length($1)}' > $TMP/${i}.readCount && \
      
      # Genome Size
      G=$(awk '{s+=$2} END {print s}' $REFERENCE.fai) && \
      
      # Number of original reads
      #N_READS=$(cat $TMP/${i}.readCount | awk '{s+=$1} END {print s}') && \
      N_READS=$(cat $TMP/${i}.readCount | wc -l) && \
      
      # Average read length
      #L=$(cat $TMP/${i}.readCount | awk '{s+=$1;m+=$1*$2} END {printf "%.0f\n", m/s}') && \
      L=$(cat $TMP/${i}.readCount | datamash median 1) && \
      L2=$(echo "scale=0; $L * 2"|bc) && \
      
      # Number of sampling reads
      READ_COUNT_DOWNSAMPLING=$(echo "scale=0; $DOWNSAMPLING_FASTQ * $G / $L2" |bc) && \
      seed=$(echo $RANDOM) && \
      
      if [[ "$N_READS" -gt "$READ_COUNT_DOWNSAMPLING" ]]; then
         if echo $EXTENSION_PE1 | grep -q "gz"; then
           NEW_PE1=$(echo $EXTENSION_PE1 | rev | cut -c 4- | rev)
           NEW_PE2=$(echo $EXTENSION_PE2 | rev | cut -c 4- | rev)
           wait
           seqtk sample -s$seed $INPUT_DIR/"${i}"$EXTENSION_PE1 $READ_COUNT_DOWNSAMPLING > $DOWNSAMPLING_FOLDER/"$i"$NEW_PE1
           seqtk sample -s$seed $INPUT_DIR/"${i}"$EXTENSION_PE2 $READ_COUNT_DOWNSAMPLING > $DOWNSAMPLING_FOLDER/"$i"$NEW_PE2
           wait
           gzip -f $DOWNSAMPLING_FOLDER/"$i"$NEW_PE1
           gzip -f $DOWNSAMPLING_FOLDER/"$i"$NEW_PE2
           wait
         else
           seqtk sample -s$seed $INPUT_DIR/"${i}"$EXTENSION_PE1 $READ_COUNT_DOWNSAMPLING > $DOWNSAMPLING_FOLDER/"$i"$EXTENSION_PE1
           seqtk sample -s$seed $INPUT_DIR/"${i}"$EXTENSION_PE2 $READ_COUNT_DOWNSAMPLING > $DOWNSAMPLING_FOLDER/"$i"$EXTENSION_PE2
         fi
      fi && \
      wait
      
      unset READ_COUNT && \
      unset READ_LENGTH && \
      unset GENOME_SIZE && \
      unset READ_COUNT_DOWNSAMPLING && \
      unset seed && \
      rm -r $TMP/${i}.readCount 2> /dev/null &
      wait
      
      # limit jobs
      if (( $(($((++n)) % $BATCH)) == 0 )) ; then
      wait # wait until all have finished (not optimal, but most times good enough)
      echo $n Files completed
      fi
      
   done
   wait

fi
wait


END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ DOWNSAMPLING FASTQ FILES TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@#############################################################"
