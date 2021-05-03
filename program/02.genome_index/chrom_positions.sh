#!/bin/bash


chunks=$BP_BY_CHROM
cat $CHROM_SIZE | awk '{print $1"\t"$2}' > ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_positions.txt.tmp1
sed -i 's/\t/\\/g' ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_positions.txt.tmp1
dos2unix ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_positions.txt.tmp1 2> /dev/null

rm ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_positions.txt 2> /dev/null
for i in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_positions.txt.tmp1); do
   chrom=$( echo $i | awk -F\\ '{print $1}')
   chromLength=$( echo $i |awk -F\\ '{print $2}')
   division=$((chromLength / chunks))
   reimander=$(expr $chromLength % $chunks)

   if [[ "$reimander" == 0 ]]; then
      chrom_rep=$(seq $division | sed "c $chrom")
      start=$(echo $(seq 1 $chunks $chromLength))
      end=$(echo $(seq $chunks $chunks $chromLength))

      paste -d ":"  <(echo "$chrom_rep" | tr -s ' '  '\n') <(echo "$start" | tr -s ' '  '\n') > ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_positions.txt.tmp2
      paste -d "-" <(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_positions.txt.tmp2) <(echo "$end" | tr -s ' '  '\n') >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_positions.txt
      rm ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_positions.txt.tmp1 2> /dev/null
      rm ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_positions.txt.tmp2 2> /dev/null
   
   else
   
      chrom_rep=$(seq $(($division+1)) | sed "c $chrom")
      start=$(echo $(seq 1 $chunks $chromLength))
      end=$(echo $(seq $chunks $chunks $chromLength) $chromLength)   

      paste -d ":"  <(echo "$chrom_rep" | tr -s ' '  '\n') <(echo "$start" | tr -s ' '  '\n') > ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_positions.txt.tmp2
      paste -d "-" <(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_positions.txt.tmp2) <(echo "$end" | tr -s ' '  '\n') >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_positions.txt
      rm ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_positions.txt.tmp1 2> /dev/null
      rm ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/tmp_${ANALYSIS_ID}/chromosomes_${ANALYSIS_ID}/chromosome_positions.txt.tmp2 2> /dev/null
    
   fi
done