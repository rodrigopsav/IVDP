#!/bin/bash


############################
##### CHROMOSOME NAMES #####
############################

echo
echo
echo "#@#############################################################"
echo "#@ GETTING CHROMOSOME NAME AND LENGTH FROM REFERENCE GENOME"
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo


### INDEX REFERENCE FILE
if test -f "$(echo $REFERENCE.fai)"; then 
   echo "INDEXED REFERENCE.FAI ALREADY EXISTS"
else
   samtools faidx $REFERENCE
   echo "INDEXED REFERENCE.FAI WAS CREATED SUCCESSFULLY"
fi
wait

### GET CHROMOSOME NAME
TOTAL_CHROM=$(echo $SELECT_CHROM | tr " " "\n" | wc -l)
if [[ "$INCLUDE_CHROM_X" == "yes" ]]; then TOTAL_CHROM=$((TOTAL_CHROM +1)); fi
if [[ "$INCLUDE_CHROM_Y" == "yes" ]]; then TOTAL_CHROM=$((TOTAL_CHROM +1)); fi
if [[ "$INCLUDE_CHROM_MT" == "yes" ]]; then TOTAL_CHROM=$((TOTAL_CHROM +1)); fi


if [[ -f ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_name.txt ]]; then

   if [[ "$TOTAL_CHROM" == "$(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_name.txt | wc -l)" ]]; then
      export CHROM_NAME=${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_name.txt
      echo "chromosome_name.txt WAS CREATED IN ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}"
   
   else
   
      rm ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_name.txt 2> /dev/null
      
      for chr in $SELECT_CHROM; do 
         sed -n ${chr}p $(echo $REFERENCE.fai) | cut -f1 >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_name.txt
      done
      wait

      if [[ "$INCLUDE_CHROM_X" == "yes" ]]; then
         sed -n $((CHROM_SET + 1))p $(echo $REFERENCE.fai) | cut -f1 >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_name.txt
      fi
      wait
      
      if [[ "$INCLUDE_CHROM_Y" == "yes" ]]; then
         sed -n $((CHROM_SET + 2))p $(echo $REFERENCE.fai) | cut -f1 >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_name.txt
      fi
      wait
      
      if [[ "$INCLUDE_CHROM_MT" == "yes" ]]; then
         sed -n $((CHROM_SET + 3))p $(echo $REFERENCE.fai) | cut -f1 >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_name.txt
      fi
      wait
      
      export CHROM_NAME=${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_name.txt
      echo "chromosome_name.txt WAS CREATED IN ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}"
   fi
   
else

   mkdir -p ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID} > /dev/null 2>&1
   for chr in $SELECT_CHROM; do 
      sed -n ${chr}p $(echo $REFERENCE.fai) | cut -f1 >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_name.txt
   done
   wait

   if [[ "$INCLUDE_CHROM_X" == "yes" ]]; then
      sed -n $((CHROM_SET + 1))p $(echo $REFERENCE.fai) | cut -f1 >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_name.txt
   fi
   wait
   
   if [[ "$INCLUDE_CHROM_Y" == "yes" ]]; then
      sed -n $((CHROM_SET + 2))p $(echo $REFERENCE.fai) | cut -f1 >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_name.txt
   fi
   wait
   
   if [[ "$INCLUDE_CHROM_MT" == "yes" ]]; then
      sed -n $((CHROM_SET + 3))p $(echo $REFERENCE.fai) | cut -f1 >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_name.txt
   fi
   wait
      
   export CHROM_NAME=${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_name.txt
   echo "chromosome_name.txt WAS CREATED IN ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}"
fi
wait


### GET CHROMOSOME LENGTH
if [[ -f ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_size.txt ]]; then

   if [[ "TOTAL_CHROM" == "${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_size.txt | wc -l)" ]]; then
      export CHROM_SIZE=${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_size.txt
      echo "chromosome_size.txt WAS CREATED IN ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}"
   
   else
   
      rm ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_size.txt 2> /dev/null
      
      for chr in $SELECT_CHROM; do 
         sed -n ${chr}p $(echo $REFERENCE.fai) | cut -f1,2 >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_size.txt
      done
      wait

      if [[ "$INCLUDE_CHROM_X" == "yes" ]]; then
         sed -n $((CHROM_SET + 1))p $(echo $REFERENCE.fai) | cut -f1,2 >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_size.txt
      fi
      wait
      
      if [[ "$INCLUDE_CHROM_Y" == "yes" ]]; then
         sed -n $((CHROM_SET + 2))p $(echo $REFERENCE.fai) | cut -f1,2 >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_size.txt
      fi
      wait
      
      if [[ "$INCLUDE_CHROM_MT" == "yes" ]]; then
         sed -n $((CHROM_SET + 3))p $(echo $REFERENCE.fai) | cut -f1,2 >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_size.txt
      fi
      wait
      
      export CHROM_SIZE=${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_size.txt
      echo "chromosome_size.txt WAS CREATED IN ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}"
   fi
   
else

   for chr in $SELECT_CHROM; do 
      sed -n ${chr}p $(echo $REFERENCE.fai) | cut -f1,2 >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_size.txt
   done
   wait

   if [[ "$INCLUDE_CHROM_X" == "yes" ]]; then
      sed -n $((CHROM_SET + 1))p $(echo $REFERENCE.fai) | cut -f1,2 >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_size.txt
   fi
   wait
   
   if [[ "$INCLUDE_CHROM_Y" == "yes" ]]; then
      sed -n $((CHROM_SET + 2))p $(echo $REFERENCE.fai) | cut -f1,2 >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_size.txt
   fi
   wait
   
   if [[ "$INCLUDE_CHROM_MT" == "yes" ]]; then
      sed -n $((CHROM_SET + 3))p $(echo $REFERENCE.fai) | cut -f1,2 >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_size.txt
   fi
   wait
   
   export CHROM_SIZE=${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_size.txt
   echo "chromosome_size.txt WAS CREATED IN ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}"
fi
wait


### SPLIT CHROMOSOME IN CHUNKS OF BASE PAIRS (BP)
rm -r ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_positions.txt 2> /dev/null
if [[ "$BP_BY_CHROM" != "all" ]]; then
   $IVDP_DIR/program/02.genome_index/chrom_positions.sh 2> /dev/null
   wait
   export CHROM_POSITIONS=${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_positions.txt
   echo "chromosome_positions.txt WAS CREATED IN ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}"
fi
wait
echo "#@#############################################################"
echo