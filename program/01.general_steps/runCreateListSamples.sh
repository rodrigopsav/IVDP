#!/bin/bash

#echo "#@#############################################################"
echo "#@ LIST OF SAMPLES"
echo

mkdir -p ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report 2>/dev/null
rm -r ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt 2>/dev/null
rm -r ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_RGID_RGPU.txt 2>/dev/null
rm -r ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt 2>/dev/null

if [[ "${LIST_SAMPLES,,}" == "none" ]]; then

   if [[ "${TYPE,,}" == "se" ]]; then
      echo "LIST_SAMPLES: Create list of samples with single-end files"
      
      
      ### SAMPLE NAMES
      cd $INPUT_DIR
      
      if [[ "$(ls *$EXTENSION_SE | wc -l)" == 0 ]]; then
         echo "Error: Check EXTENSION_SE or INPUT_DIR variables in the parameter card"
         echo "Aborting analysis"
         exit 1
      fi
      
      ls *$EXTENSION_SE >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt
      sed -i "s/$EXTENSION_SE//g" ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt

      ### RGID and RGPU:
      for SAMPLE in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt); do
         HEADER=$(zcat -f ${INPUT_DIR}/${SAMPLE}${EXTENSION_SE} | head -1)
         
         if echo "$HEADER" | grep -q "SRR\|ERR"; then
            echo $SAMPLE >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_RGID_RGPU.txt
         else
            echo $HEADER | awk -F ":" '{ print $3"."$4"."$10 }' >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_RGID_RGPU.txt
         fi
      done

      ### JOIT THE FILES
      paste -d "\t" ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt \
                    ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt \
                    ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_RGID_RGPU.txt \
                    >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt
      sed -i $'1 i\\\nSAMPLE\tINDIVIDUAL_ID\tREAD_GROUP_ID' ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt

      rm -r ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt 2>/dev/null
      rm -r ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_RGID_RGPU.txt 2>/dev/null


   else
      echo "LIST_SAMPLES: Create list of samples with paired-end files"


      ### SAMPLE NAMES
      cd $INPUT_DIR
      
      if [[ "$(ls *$EXTENSION_PE1 | wc -l)" == 0 || "$(ls *$EXTENSION_PE2 | wc -l)" == 0 ]]; then
         echo "Error: Check EXTENSION_PE1, EXTENSION_PE2 or INPUT_DIR variables in the parameter card"
         echo "Aborting analysis"
         exit 1
      fi
      
      ls *$EXTENSION_PE1 >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt
      sed -i "s/$EXTENSION_PE1//g" ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt

      ### RGID and RGPU:
      for SAMPLE in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt); do
         HEADER=$(zcat -f ${INPUT_DIR}/${SAMPLE}${EXTENSION_PE1} | head -1)
         
         if echo "$HEADER" | grep -q "SRR\|ERR"; then
            echo $SAMPLE >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_RGID_RGPU.txt
         else
            echo $HEADER | awk -F ":" '{ print $3"."$4"."$10 }' >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_RGID_RGPU.txt
         fi
      done

      ### JOIT THE FILES
      paste -d "\t" ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt \
                    ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt \
                    ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_RGID_RGPU.txt \
                    >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt
      sed -i $'1 i\\\nSAMPLE\tINDIVIDUAL_ID\tREAD_GROUP_ID' ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt

      rm -r ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt 2>/dev/null
      rm -r ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_RGID_RGPU.txt 2>/dev/null

   fi

else
   
   if [[ "${TYPE,,}" == "se" ]]; then
      echo "LIST_SAMPLES: Create list of samples with single-end files"
      echo "LIST_SAMPLES in $(readlink -f $LIST_SAMPLES)"      
      
      ### CHANGE ONE OR MULTIPLE SPACES BY TAB
      # https://serverfault.com/questions/431167/sed-replace-all-tabs-and-spaces-with-a-single-space/431168
      sed -e "s/[[:space:]]\+/\t/g" $LIST_SAMPLES >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt
      sed -i 's/,/\t/g' ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt
      awk -F'\t' '{print $1"\t"$2}' ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt > ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt.bak
      mv ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt.bak ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt
      sed -i 's/\t/\\/g' ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt
      
      ### RGID and RGPU:
      for line in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt); do
          SAMPLE=$(echo $line | awk -F\\ '{ print $1}')
          INDIVIDUAL_ID=$(echo $line | awk -F\\ '{print $2}')
          HEADER=$(zcat -f ${INPUT_DIR}/${SAMPLE}${EXTENSION_SE} | head -1)
          
          if echo "$HEADER" | grep -q "SRR\|ERR"; then
             echo "$SAMPLE"_"$INDIVIDUAL_ID" >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_RGID_RGPU.txt
          else
             echo $HEADER | awk -F ":" '{ print $3"."$4"."$10 }' >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_RGID_RGPU.txt
          fi
      done

      ### JOIN THE FILES
      paste -d "\t" ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt \
                    ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_RGID_RGPU.txt \
                    >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt
      sed -i $'1 i\\\nSAMPLE\tINDIVIDUAL_ID\tREAD_GROUP_ID' ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt
      
      rm -r ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt 2>/dev/null
      rm -r ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_RGID_RGPU.txt 2>/dev/null
      sed -i 's/\\/\t/g' ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt
      

   else
      echo "LIST_SAMPLES: Create list of samples with paired-end files"
      echo "LIST_SAMPLES in $(readlink -f $LIST_SAMPLES)" 

      ### CHANGE ONE OR MULTIPLE SPACES BY TAB
      sed -e "s/[[:space:]]\+/\t/g" $LIST_SAMPLES >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt
      sed -i 's/,/\t/g' ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt
      awk -F'\t' '{print $1"\t"$2}' ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt > ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt.bak
      mv ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt.bak ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt
      sed -i 's/\t/\\/g' ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt
      
      ### RGID and RGPU:
      for line in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt); do
          SAMPLE=$(echo $line | awk -F\\ '{ print $1}')
          INDIVIDUAL_ID=$(echo $line | awk -F\\ '{print $2}')
          HEADER=$(zcat -f ${INPUT_DIR}/${SAMPLE}${EXTENSION_PE1} | head -1)
          
          if echo "$HEADER" | grep -q "SRR\|ERR"; then
             echo "$SAMPLE"_"$INDIVIDUAL_ID" >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_RGID_RGPU.txt
          else
             echo $HEADER | awk -F ":" '{ print $3"."$4"."$10 }' >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_RGID_RGPU.txt
          fi
      done

      ### JOIN THE FILES
      paste -d "\t" ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt \
                    ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_RGID_RGPU.txt \
                    >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt
      sed -i $'1 i\\\nSAMPLE\tINDIVIDUAL_ID\tREAD_GROUP_ID' ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt
      
      rm -r ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_samples.txt 2>/dev/null
      rm -r ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_list_RGID_RGPU.txt 2>/dev/null
      sed -i 's/\\/\t/g' ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt
      
   fi
fi

#echo
echo "#@#############################################################"

