#!/bin/bash

echo
echo
echo "#@#############################################################"
echo "#@                 GENE_COUNT: $OUTPUT_NAME "
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo

START1=$(date +%s)

######################
##### GENE COUNT #####
######################


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


##### Check Output Folders for This Step
for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
   for GENE_COUNT in $GENECOUNT_PROGRAM; do
      export GENE_COUNT=$GENE_COUNT
      for FEATURE in $FEATURE_TYPE; do
         export FEATURE=$FEATURE
      
         source "$IVDP_DIR"/program/01.general_steps/folders.sh
   
         mkdir -p "$GENECOUNT_DATA" > /dev/null 2>&1
         mkdir -p "$LOG_GENECOUNT" > /dev/null 2>&1
         mkdir -p "$REPORT_GENECOUNT" > /dev/null 2>&1
        
      done
   done      
done
wait


##### Run Gene Count for RNAseq samples
n=0
for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
   for GENE_COUNT in $GENECOUNT_PROGRAM; do
      export GENE_COUNT=$GENE_COUNT
      for FEATURE in $FEATURE_TYPE; do
         export FEATURE=$FEATURE
   	
         for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
            export i=$( echo $sample | awk -F\\ '{print $1}')
            export SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
            export ID_PU=$( echo $sample |awk -F\\ '{print $3}')
            echo "$ALIGNER $GENE_COUNT $FEATURE: SAMPLE $i"
            source "$IVDP_DIR"/program/01.general_steps/folders.sh
   
   
            ##### GENE COUNTS
		    $IVDP_DIR/program/05.gene_count/genecount_"$GENE_COUNT".sh >> $LOG_GENECOUNT/"$i".txt 2>&1 &
   
   
		    # limit jobs
		    if (( $(($((++n)) % $BATCH)) == 0 )) ; then
		    wait # wait until all have finished (not optimal, but most times good enough)
		    echo $n Files completed
		    fi
         done
      done
   done
done
wait
echo


##### Create Gene Count Tables with All the Samples
for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
   for GENE_COUNT in $GENECOUNT_PROGRAM; do
      export GENE_COUNT=$GENE_COUNT
      for FEATURE in $FEATURE_TYPE; do
         export FEATURE=$FEATURE
         source "$IVDP_DIR"/program/01.general_steps/folders.sh         
         
         cd $GENECOUNT_DATA
         if [[ "$GENE_COUNT" == "htseq" ]]; then
            echo -e "Geneid\t$(ls *$FEATURE.txt | rev | cut -d_ -f3- | rev | tr "\n" "\t" | sed 's/.$//')" > sampleList.txt
            cat $(ls *$FEATURE.txt | head -1) | head -n -5 | awk '{print $1}' > ensemblId.txt
            cat ensemblId.txt > "$GENE_COUNT"_"$FEATURE"_"$ALIGNER"_CountTable.txt
            
            for file in $( ls *$FEATURE.txt ); do
               cat $file | head -n -5 > tmp.txt
               awk '{print $2}' < tmp.txt | paste "$GENE_COUNT"_"$FEATURE"_"$ALIGNER"_CountTable.txt - >> "$GENE_COUNT"_"$FEATURE"_"$ALIGNER"_CountTable.txt.tmp && mv "$GENE_COUNT"_"$FEATURE"_"$ALIGNER"_CountTable.txt.tmp "$GENE_COUNT"_"$FEATURE"_"$ALIGNER"_CountTable.txt
            done         
            
            cat sampleList.txt "$GENE_COUNT"_"$FEATURE"_"$ALIGNER"_CountTable.txt > "$GENE_COUNT"_"$FEATURE"_"$ALIGNER"_CountTable.txt.tmp && mv "$GENE_COUNT"_"$FEATURE"_"$ALIGNER"_CountTable.txt.tmp "$GENE_COUNT"_"$FEATURE"_"$ALIGNER"_CountTable.txt
            
            rm sampleList.txt
            rm ensemblId.txt
            rm tmp.txt
            mv "$GENE_COUNT"_"$FEATURE"_"$ALIGNER"_CountTable.txt ${GENECOUNT_DATA%/*}
            
         else
            
            echo -e "Geneid\t$(ls *$FEATURE.txt | rev | cut -d_ -f3- | rev | tr "\n" "\t" | sed 's/.$//')" > sampleList.txt
            grep -v "#" $(ls *$FEATURE.txt | head -1) | tail -n +2 | awk '{print $1}' > ensemblId.txt
            cat ensemblId.txt > "$GENE_COUNT"_"$FEATURE"_"$ALIGNER"_CountTable.txt
            
            for file in $( ls *$FEATURE.txt ); do
               grep -v "#" $file | tail -n +2 > tmp.txt
               awk '{print $7}' < tmp.txt | paste "$GENE_COUNT"_"$FEATURE"_"$ALIGNER"_CountTable.txt - >> "$GENE_COUNT"_"$FEATURE"_"$ALIGNER"_CountTable.txt.tmp && mv "$GENE_COUNT"_"$FEATURE"_"$ALIGNER"_CountTable.txt.tmp "$GENE_COUNT"_"$FEATURE"_"$ALIGNER"_CountTable.txt
            done
            
            cat sampleList.txt "$GENE_COUNT"_"$FEATURE"_"$ALIGNER"_CountTable.txt > "$GENE_COUNT"_"$FEATURE"_"$ALIGNER"_CountTable.txt.tmp && mv "$GENE_COUNT"_"$FEATURE"_"$ALIGNER"_CountTable.txt.tmp "$GENE_COUNT"_"$FEATURE"_"$ALIGNER"_CountTable.txt

            rm sampleList.txt
            rm ensemblId.txt
            rm tmp.txt
            mv "$GENE_COUNT"_"$FEATURE"_"$ALIGNER"_CountTable.txt ${GENECOUNT_DATA%/*}
            
         fi
         
      done
   done
done
wait


##############################
##### MULTIQC GENE_COUNT #####
##############################
echo "MULTIQC REPORT"
	
for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
   for GENE_COUNT in $GENECOUNT_PROGRAM; do
      export GENE_COUNT=$GENE_COUNT
      for FEATURE in $FEATURE_TYPE; do
         export FEATURE=$FEATURE
         source "$IVDP_DIR"/program/01.general_steps/folders.sh
      
         $IVDP_DIR/program/05.gene_count/genecount_multiqc.sh > /dev/null 2>&1 &
      
      done
   done
done
wait

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ GENE COUNT TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@#############################################################"
