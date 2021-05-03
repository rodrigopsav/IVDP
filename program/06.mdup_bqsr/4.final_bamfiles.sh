#!/bin/bash


# OPERATORS IN BASH
# https://www.geeksforgeeks.org/basic-operators-in-shell-scripting/

# IF STATEMENT
# https://ryanstutorials.net/bash-scripting-tutorial/bash-if-statements.php

# EXPRESSIONS IN BASH IN IF STATEMENT
# https://www.tutorialkart.com/bash-shell-scripting/bash-if/

# IF STATEMENT AND COMPARISON OPERATORS
# https://www.sanspire.com/bash-if-statement-and-comparison-operators/

# TEST CONDITIONS IN BASH
# https://www.lifewire.com/test-linux-command-unix-command-4097166


splitncigar="$(ls $ALIGN_MDUP/*splitncigar.bam 2> /dev/null | wc -l)"
echo "#@ NUMBER OF SPLITNCIGAR FILES = $splitncigar"

mdup="$(ls $ALIGN_MDUP/*RG_mdup.bam 2> /dev/null | wc -l)"
echo "#@ NUMBER OF MDUP FILES = $mdup"

bqsr="$(ls $ALIGN_BQSR/*bqsr.bam 2> /dev/null | wc -l)"
echo "#@ NUMBER OF BQSR FILES = $bqsr"


##### MDUP AND BQSR FILES
if [[ "$splitncigar" != 0 ]]; then
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      export i=$( echo $sample | awk -F\\ '{print $1}')
      export SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      export ID_PU=$( echo $sample |awk -F\\ '{print $3}')
      
      rm -r $ALIGN_MDUP/${i}_trim_"$ALIGNER"_sorted_RG_mdup.bam 2> /dev/null
      rm -r $ALIGN_MDUP/${i}_trim_"$ALIGNER"_sorted_RG_mdup.bam.bai 2> /dev/null
      
      mv $ALIGN_MDUP/${i}_trim_"$ALIGNER"_sorted_RG_mdup_splitncigar.bam $ALIGN_MDUP/${i}_trim_"$ALIGNER"_sorted_RG_mdup.bam 2> /dev/null
      mv $ALIGN_MDUP/${i}_trim_"$ALIGNER"_sorted_RG_mdup_splitncigar.bai $ALIGN_MDUP/${i}_trim_"$ALIGNER"_sorted_RG_mdup.bai 2> /dev/null
   done
fi
wait


if [[ "$mdup" != 0 ]] && [[ "$bqsr" != 0 ]]; then
   echo "#@ MDUP, AND BQSR files exist"
   echo "#@ Use $ALIGN_BQSR/*bqsr.bam for variant calling (except for FREEBAYES)"
   echo "#@ Use $ALIGN_MDUP/*mdup.bam for variant calling with FREEBAYES"


##### RG FILES
elif [[ "$mdup" != 0 ]] && [[ "$bqsr" == 0 ]]; then
   echo "#@ ONLY MDUP files exist"
   echo "#@ Use $ALIGN_MDUP/*mdup.bam for variant calling"

else
   echo "#@ ERROR: THERE ARE NO MDUP OR BQSR BAM FILES"
   echo "#@ PLEASE CHECK IF THE BAM FILES FROM ALIGNMENT STEP ARE OK (*RG.bam files in $ALIGN SUBFOLDERS)"
   echo "#@ PLEASE CHECK THE LOG FILES IN:"
   echo "$LOG_ALIGN"
   echo "$LOG_ALIGN_MDUP"
   echo "$LOG_ALIGN_BQSR"
   exit 1
fi

