#!/bin/bash

#IMPORTANT: NEVER LEAVE A SPACE AFTER ANY \ WHEN USING SBATCH

#-------------------------------#
##### STEP1: ASSEMBLY INDEX #####
#-------------------------------#
rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid1.txt > /dev/null 2>&1
JOBNAME1="${ANALYSIS_ID}"_01.Index

JOBID1=$(sbatch --job-name=$JOBNAME1 --time=08:00:00 --cpus-per-task=1 --mem=200G \
--output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME1".txt \
--export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME \
$IVDP_DIR/program/01.general_steps/slurm/1a.slurmAssemblyidx.sb)
echo $JOBID1 >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid1.txt
wait

WAIT1=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid1.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
wait


#-------------------------------#
##### STEP1b: FASTQ SUMMARY #####
#-------------------------------#
rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid1b.txt > /dev/null 2>&1
JOBNAME1b="${ANALYSIS_ID}"_01b.FastqSummary

JOBID1b=$(sbatch --job-name=$JOBNAME1b --time=24:00:00 --cpus-per-task=32 --mem=8G \
--output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME1b".txt \
--export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME \
$IVDP_DIR/program/01.general_steps/slurm/1b.slurmFastqStats.sb)
echo $JOBID1b >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid1b.txt
wait

WAIT1b=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid1b.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
wait


#-----------------------------------#
##### STEP2: DOWNSAMPLING FASTQ #####
#-----------------------------------#
if [[ $DOWNSAMPLING_FASTQ != 0 ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid2.txt > /dev/null 2>&1
   
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      i=$( echo $sample | awk -F\\ '{print $1}')
      SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      ID_PU=$( echo $sample |awk -F\\ '{print $3}')
      JOBNAME2="${ANALYSIS_ID}"_02.DSFastq_"${i}"
      
      JOBID2=$(sbatch --dependency=afterok:$WAIT1 --job-name=$JOBNAME2 --time=08:00:00 --cpus-per-task=1 --mem=80G \
      --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME2".txt \
      --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,i=$i,SUBJECT_ID=$SUBJECT_ID,ID_PU=$ID_PU \
      $IVDP_DIR/program/01.general_steps/slurm/2.slurmDownSamplingFastq.sb)
      echo $JOBID2 >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid2.txt
   done
   wait
   
   WAIT2=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid2.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
fi
wait


#-----------------------------------#
###### STEP3: ADAPTER TRIMMING ######
#-----------------------------------#
if [[ $DOWNSAMPLING_FASTQ != 0 && $STEP_QC_TRIM == "yes" ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid3.txt > /dev/null 2>&1
   
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      i=$( echo $sample | awk -F\\ '{print $1}')
      SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      ID_PU=$( echo $sample |awk -F\\ '{print $3}')
      JOBNAME3="${ANALYSIS_ID}"_03.Trim_"${i}"
   
      JOBID3=$(sbatch --dependency=afterok:$WAIT1:$WAIT2 --job-name=$JOBNAME3 --time=$TRIMMING_TIME --cpus-per-task=$TRIMMING_CPU --mem=$TRIMMING_MEM \
      --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME3".txt \
      --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,i=$i,SUBJECT_ID=$SUBJECT_ID,ID_PU=$ID_PU \
      $IVDP_DIR/program/01.general_steps/slurm/3.slurmTrim.sb)
      echo $JOBID3 >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid3.txt
   done
   wait
   
   WAIT3=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid3.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
fi
wait


if [[ $DOWNSAMPLING_FASTQ == 0 && $STEP_QC_TRIM == "yes" ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid3.txt > /dev/null 2>&1
   
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      i=$( echo $sample | awk -F\\ '{print $1}')
      SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      ID_PU=$( echo $sample |awk -F\\ '{print $3}')
      JOBNAME3="${ANALYSIS_ID}"_03.Trim_"${i}"
   
      JOBID3=$(sbatch --dependency=afterok:$WAIT1 --job-name=$JOBNAME3 --time=$TRIMMING_TIME --cpus-per-task=$TRIMMING_CPU --mem=$TRIMMING_MEM \
      --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME3".txt \
      --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,i=$i,SUBJECT_ID=$SUBJECT_ID,ID_PU=$ID_PU \
      $IVDP_DIR/program/01.general_steps/slurm/3.slurmTrim.sb)
      echo $JOBID3 >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid3.txt
   done
   wait
   
   WAIT3=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid3.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
fi
wait


#----------------------------#
###### STEP4: ALIGNMENT ######
#----------------------------#
if [[ $DOWNSAMPLING_FASTQ != 0 && $STEP_QC_TRIM == "yes" && $STEP_ALIGNMENT == "yes" ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid4.txt > /dev/null 2>&1
   
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      i=$( echo $sample | awk -F\\ '{print $1}')
      SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      ID_PU=$( echo $sample |awk -F\\ '{print $3}')
      JOBNAME4="${ANALYSIS_ID}"_04.Align_"${i}"
         
      JOBID4=$(sbatch --dependency=afterok:$WAIT1:$WAIT2:$WAIT3 --job-name=$JOBNAME4 --time=$ALIGNMENT_TIME --cpus-per-task=$ALIGNMENT_CPU --mem=$ALIGNMENT_MEM \
      --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME4".txt \
      --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,i=$i,SUBJECT_ID=$SUBJECT_ID,ID_PU=$ID_PU \
      $IVDP_DIR/program/01.general_steps/slurm/4a.slurmAlignment.sb)
      echo $JOBID4 >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid4.txt
   done
   wait
   
   WAIT4=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid4.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
fi
wait


if [[ $DOWNSAMPLING_FASTQ == 0 && $STEP_QC_TRIM == "yes" && $STEP_ALIGNMENT == "yes" ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid4.txt > /dev/null 2>&1
   
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      i=$( echo $sample | awk -F\\ '{print $1}')
      SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      ID_PU=$( echo $sample |awk -F\\ '{print $3}')
      JOBNAME4="${ANALYSIS_ID}"_04.Align_"${i}"
         
      JOBID4=$(sbatch --dependency=afterok:$WAIT1:$WAIT3 --job-name=$JOBNAME4 --time=$ALIGNMENT_TIME --cpus-per-task=$ALIGNMENT_CPU --mem=$ALIGNMENT_MEM \
      --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME4".txt \
      --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,i=$i,SUBJECT_ID=$SUBJECT_ID,ID_PU=$ID_PU \
      $IVDP_DIR/program/01.general_steps/slurm/4a.slurmAlignment.sb)
      echo $JOBID4 >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid4.txt
   done
   wait
   
   WAIT4=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid4.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
fi
wait


if [[ $DOWNSAMPLING_FASTQ != 0 && $STEP_QC_TRIM == "no" && $STEP_ALIGNMENT == "yes" ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid4.txt > /dev/null 2>&1
   
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      i=$( echo $sample | awk -F\\ '{print $1}')
      SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      ID_PU=$( echo $sample |awk -F\\ '{print $3}')
      JOBNAME4="${ANALYSIS_ID}"_04.Align_"${i}"
         
      JOBID4=$(sbatch --dependency=afterok:$WAIT1:$WAIT2 --job-name=$JOBNAME4 --time=$ALIGNMENT_TIME --cpus-per-task=$ALIGNMENT_CPU --mem=$ALIGNMENT_MEM \
      --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME4".txt \
      --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,i=$i,SUBJECT_ID=$SUBJECT_ID,ID_PU=$ID_PU \
      $IVDP_DIR/program/01.general_steps/slurm/4a.slurmAlignment.sb)
      echo $JOBID4 >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid4.txt
   done
   wait
   
   WAIT4=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid4.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
fi
wait


if [[ $DOWNSAMPLING_FASTQ == 0 && $STEP_QC_TRIM == "no" && $STEP_ALIGNMENT == "yes" ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid4.txt > /dev/null 2>&1
   
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      i=$( echo $sample | awk -F\\ '{print $1}')
      SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      ID_PU=$( echo $sample |awk -F\\ '{print $3}')
      JOBNAME4="${ANALYSIS_ID}"_04.Align_"${i}"
         
      JOBID4=$(sbatch --dependency=afterok:$WAIT1 --job-name=$JOBNAME4 --time=$ALIGNMENT_TIME --cpus-per-task=$ALIGNMENT_CPU --mem=$ALIGNMENT_MEM \
      --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME4".txt \
      --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,i=$i,SUBJECT_ID=$SUBJECT_ID,ID_PU=$ID_PU \
      $IVDP_DIR/program/01.general_steps/slurm/4a.slurmAlignment.sb)
      echo $JOBID4 >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid4.txt
   done
   wait
   
   WAIT4=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid4.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
fi
wait


#----------------------------#
##### STEP4b: STATISTICS #####
#----------------------------#
if [[ $STEP_ALIGNMENT == "yes" && $STATS_ALIGNMENT == "yes" ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid4b.txt > /dev/null 2>&1
   JOBNAME4b="${ANALYSIS_ID}"_04b.Stats_Align
   
   JOBID4b=$(sbatch --dependency=afterok:$WAIT4 --job-name=$JOBNAME4b \
   --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME4b".txt \
   --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME \
   $IVDP_DIR/program/01.general_steps/slurm/4b.slurmStatsAlignment.sb)
   echo $JOBID4b >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid4b.txt
   
   WAIT4b=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid4b.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
fi
wait


#-----------------------------#
###### STEP5: GENE COUNT ######
#-----------------------------#
##### Gene Counts
if [[ $STEP_ALIGNMENT == "yes" && $STEP_GENECOUNT == "yes" && $GENECOUNT_PROGRAM != "none" ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid5a.txt > /dev/null 2>&1
   
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      i=$( echo $sample | awk -F\\ '{print $1}')
      SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      ID_PU=$( echo $sample |awk -F\\ '{print $3}')
      JOBNAME5a="${ANALYSIS_ID}"_05.GeneCount_"${i}"
         
      JOBID5a=$(sbatch --dependency=afterok:$WAIT1:$WAIT4 --job-name=$JOBNAME5a --time=$GENECOUNT_TIME --cpus-per-task=$GENECOUNT_CPU --mem=$GENECOUNT_MEM \
      --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME5a".txt \
      --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,i=$i,SUBJECT_ID=$SUBJECT_ID,ID_PU=$ID_PU \
      $IVDP_DIR/program/01.general_steps/slurm/5a.slurmGeneCount.sb)
      echo $JOBID5a >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid5a.txt
   done
   wait
   
   WAIT5a=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid5a.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
fi
wait


if [[ $STEP_ALIGNMENT == "no" && $STEP_GENECOUNT == "yes" && $GENECOUNT_PROGRAM != "none" ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid5a.txt > /dev/null 2>&1
   
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      i=$( echo $sample | awk -F\\ '{print $1}')
      SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      ID_PU=$( echo $sample |awk -F\\ '{print $3}')
      JOBNAME5a="${ANALYSIS_ID}"_05.GeneCount_"${i}"
      
      JOBID5a=$(sbatch --dependency=afterok:$WAIT1 --job-name=$JOBNAME5a --time=$GENECOUNT_TIME --cpus-per-task=$GENECOUNT_CPU --mem=$GENECOUNT_MEM \
      --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME5a".txt \
      --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,i=$i,SUBJECT_ID=$SUBJECT_ID,ID_PU=$ID_PU \
      $IVDP_DIR/program/01.general_steps/slurm/5a.slurmGeneCount.sb)
      echo $JOBID5a >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid5a.txt
   done
   wait
   
   WAIT5a=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid5a.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
fi
wait


##### Create Table with Gene Counts
if [[ $STEP_GENECOUNT == "yes" && $GENECOUNT_PROGRAM != "none" ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid5b.txt > /dev/null 2>&1
   JOBNAME5b="${ANALYSIS_ID}"_05.GeneCountTable
   
   JOBID5b=$(sbatch --dependency=afterok:$WAIT5a --job-name=$JOBNAME5b \
   --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME5b".txt \
   --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME \
   $IVDP_DIR/program/01.general_steps/slurm/5b.slurmGeneCountTable.sb)
   echo $JOBID5b >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid5b.txt
   wait
   
   WAIT5b=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid5b.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
   wait
   WAIT5=$(echo $WAIT5a:$WAIT5b)
   wait
fi
wait
   

#---------------------------------#
##### STEP6: DOWNSAMPLING BAM #####
#---------------------------------#
if [[ $STEP_ALIGNMENT == "yes" && $DOWNSAMPLING_BAM != 0 ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid6.txt > /dev/null 2>&1
   
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      i=$( echo $sample | awk -F\\ '{print $1}')
      SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      ID_PU=$( echo $sample |awk -F\\ '{print $3}')
      JOBNAME6="${ANALYSIS_ID}"_06.DSBam_"${i}"
         
      JOBID6=$(sbatch --dependency=afterok:$WAIT1:$WAIT4:$WAIT4b --job-name=$JOBNAME6 --time=08:00:00 --cpus-per-task=8 --mem=80G \
      --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME6".txt \
      --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,i=$i,SUBJECT_ID=$SUBJECT_ID,ID_PU=$ID_PU \
      $IVDP_DIR/program/01.general_steps/slurm/6a.slurmDownSamplingBam.sb)
      echo $JOBID6 >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid6.txt
   done
   wait
   
   WAIT6=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid6.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
fi
wait


if [[ $STEP_ALIGNMENT == "no" && $DOWNSAMPLING_BAM != 0 ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid6.txt > /dev/null 2>&1
   
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      i=$( echo $sample | awk -F\\ '{print $1}')
      SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      ID_PU=$( echo $sample |awk -F\\ '{print $3}')
      JOBNAME6="${ANALYSIS_ID}"_06.DSBam_"${i}"
         
      JOBID6=$(sbatch --dependency=afterok:$WAIT1 --job-name=$JOBNAME6 --time=$ALIGNMENT_TIME --cpus-per-task=$ALIGNMENT_CPU --mem=$ALIGNMENT_MEM \
      --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME6".txt \
      --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,i=$i,SUBJECT_ID=$SUBJECT_ID,ID_PU=$ID_PU \
      $IVDP_DIR/program/01.general_steps/slurm/6a.slurmDownSamplingBam.sb)
      echo $JOBID6 >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid6.txt
   done
   wait
   
   WAIT6=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid6.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
fi
wait


#---------------------------------------------#
##### STEP6b: STATISTICS DOWNSAMPLING BAM #####
#---------------------------------------------#
if [[ $DOWNSAMPLING_BAM != 0 ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid6b.txt > /dev/null 2>&1
   JOBNAME6b="${ANALYSIS_ID}"_6b.Stats_DSBam
   
   JOBID6b=$(sbatch --dependency=afterok:$WAIT6 --job-name=$JOBNAME6b \
   --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME6b".txt \
   --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME \
   $IVDP_DIR/program/01.general_steps/slurm/6b.slurmStatsDownSamplingBam.sb)
   echo $JOBID6b >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid6b.txt
   
   WAIT6b=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid6b.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
fi
wait


#----------------------------------------------#
###### STEP7: ALIGNMENT CORRECTION (MDUP) ######
#----------------------------------------------#
if [[ $STEP_ALIGNMENT == "yes" && $DOWNSAMPLING_BAM != 0 && $STEP_MDUP == "yes" ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid7.txt > /dev/null 2>&1
   
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      i=$( echo $sample | awk -F\\ '{print $1}')
      SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      ID_PU=$( echo $sample |awk -F\\ '{print $3}')
      JOBNAME7="${ANALYSIS_ID}"_07.Mdup_"${i}"
         
      JOBID7=$(sbatch --dependency=afterok:$WAIT1:$WAIT4:$WAIT6 --job-name=$JOBNAME7 --time=$MDUP_TIME --cpus-per-task=$MDUP_CPU --mem=$MDUP_MEM \
      --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME7".txt \
      --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,i=$i,SUBJECT_ID=$SUBJECT_ID,ID_PU=$ID_PU \
      $IVDP_DIR/program/01.general_steps/slurm/7.slurmMDUP.sb)
      echo $JOBID7 >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid7.txt
   done
   wait
   
   WAIT7=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid7.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
fi
wait


if [[ $STEP_ALIGNMENT == "no" && $DOWNSAMPLING_BAM != 0 && $STEP_MDUP == "yes" ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid7.txt > /dev/null 2>&1
   
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      i=$( echo $sample | awk -F\\ '{print $1}')
      SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      ID_PU=$( echo $sample |awk -F\\ '{print $3}')
      JOBNAME7="${ANALYSIS_ID}"_07.Mdup_"${i}"
         
      JOBID7=$(sbatch --dependency=afterok:$WAIT1:$WAIT6 --job-name=$JOBNAME7 --time=$MDUP_TIME --cpus-per-task=$MDUP_CPU --mem=$MDUP_MEM \
      --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME7".txt \
      --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,i=$i,SUBJECT_ID=$SUBJECT_ID,ID_PU=$ID_PU \
      $IVDP_DIR/program/01.general_steps/slurm/7.slurmMDUP.sb)
      echo $JOBID7 >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid7.txt
   done
   wait
   
   WAIT7=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid7.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
fi
wait


if [[ $STEP_ALIGNMENT == "yes" && $DOWNSAMPLING_BAM == 0 && $STEP_MDUP == "yes" ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid7.txt > /dev/null 2>&1
   
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      i=$( echo $sample | awk -F\\ '{print $1}')
      SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      ID_PU=$( echo $sample |awk -F\\ '{print $3}')
      JOBNAME7="${ANALYSIS_ID}"_07.Mdup_"${i}"
         
      JOBID7=$(sbatch --dependency=afterok:$WAIT1:$WAIT4 --job-name=$JOBNAME7 --time=$MDUP_TIME --cpus-per-task=$MDUP_CPU --mem=$MDUP_MEM \
      --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME7".txt \
      --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,i=$i,SUBJECT_ID=$SUBJECT_ID,ID_PU=$ID_PU \
      $IVDP_DIR/program/01.general_steps/slurm/7.slurmMDUP.sb)
      echo $JOBID7 >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid7.txt
   done
   wait
   
   WAIT7=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid7.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
fi
wait


if [[ $STEP_ALIGNMENT == "no" && $DOWNSAMPLING_BAM == 0 && $STEP_MDUP == "yes" ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid7.txt > /dev/null 2>&1
   
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      i=$( echo $sample | awk -F\\ '{print $1}')
      SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      ID_PU=$( echo $sample |awk -F\\ '{print $3}')
      JOBNAME7="${ANALYSIS_ID}"_07.Mdup_"${i}"
         
      JOBID7=$(sbatch --dependency=afterok:$WAIT1 --job-name=$JOBNAME7 --time=$MDUP_TIME --cpus-per-task=$MDUP_CPU --mem=$MDUP_MEM \
      --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME7".txt \
      --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,i=$i,SUBJECT_ID=$SUBJECT_ID,ID_PU=$ID_PU \
      $IVDP_DIR/program/01.general_steps/slurm/7.slurmMDUP.sb)
      echo $JOBID7 >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid7.txt
   done
   wait
   
   WAIT7=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid7.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
fi
wait


#----------------------------------------------#
###### STEP8: ALIGNMENT CORRECTION (BQSR) ######
#----------------------------------------------#
if [[ $STEP_MDUP == "yes" && $STEP_BQSR == "yes" ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid8.txt > /dev/null 2>&1
   
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      i=$( echo $sample | awk -F\\ '{print $1}')
      SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      ID_PU=$( echo $sample |awk -F\\ '{print $3}')
      JOBNAME8="${ANALYSIS_ID}"_08.BQSR_"${i}"
      
      JOBID8=$(sbatch --dependency=afterok:$WAIT1:$WAIT7 --job-name=$JOBNAME8 --time=$BQSR_TIME --cpus-per-task=$BQSR_CPU --mem=$BQSR_MEM \
      --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME8".txt \
      --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,i=$i,SUBJECT_ID=$SUBJECT_ID,ID_PU=$ID_PU \
      $IVDP_DIR/program/01.general_steps/slurm/8.slurmBQSR.sb)
      echo $JOBID8 >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid8.txt
   done
   wait
   
   WAIT8=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid8.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
fi
wait


if [[ $STEP_MDUP == "no" && $STEP_BQSR == "yes" ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid8.txt > /dev/null 2>&1
   
   for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
      i=$( echo $sample | awk -F\\ '{print $1}')
      SUBJECT_ID=$( echo $sample |awk -F\\ '{print $2}')
      ID_PU=$( echo $sample |awk -F\\ '{print $3}')
      JOBNAME8="${ANALYSIS_ID}"_08.BQSR_"${i}"
      
      JOBID8=$(sbatch --dependency=afterok:$WAIT1 --job-name=$JOBNAME8 --time=$BQSR_TIME --cpus-per-task=$BQSR_CPU --mem=$BQSR_MEM \
      --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME8".txt \
      --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,i=$i,SUBJECT_ID=$SUBJECT_ID,ID_PU=$ID_PU \
      $IVDP_DIR/program/01.general_steps/slurm/8.slurmBQSR.sb)
      echo $JOBID8 >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid8.txt
   done
   wait
   
   WAIT8=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid8.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
fi
wait


#---------------------------------#
##### STEP9a: VARIANT CALLING #####
#---------------------------------#
if [[ $STEP_BQSR == "yes" && $STEP_VAR_CALL == "yes" ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9a.txt > /dev/null 2>&1
   
   for CALLER in $CALLER_PROGRAM; do
      ### VCF MODE
      if [[ $CALLER != "gatk4CombineGVCFs" && $CALLER != "gatk4GenomicsDBImport" && $CALLER != "gatk3" ]]; then
         for CHR in $(cat $CHROM_NAME); do
            JOBNAME9a="${ANALYSIS_ID}"_09.CallVCFs_"${CALLER}"_"${CHR}"
         
            JOBID9a=$(sbatch --dependency=afterok:$WAIT1:$WAIT8 --job-name=$JOBNAME9a --time=$CALLVARIANT_TIME --cpus-per-task=$CALLVARIANT_CPU --mem=$CALLVARIANT_MEM \
            --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME9a".txt \
            --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,CALLER=$CALLER,chr=$CHR \
            $IVDP_DIR/program/01.general_steps/slurm/9a.slurmCallingVCFs.sb)
            echo $JOBID9a >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9a.txt
         done
         wait
      ### CALL GVCF MODE
      else  
         for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
            i=$( echo $sample | awk -F\\ '{print $1}')
            
            #for CHR in $(cat $CHROM_NAME); do
            JOBNAME9a="${ANALYSIS_ID}"_09.CallGVCFs_"${CALLER}"_"${i}"
            
            JOBID9a=$(sbatch --dependency=afterok:$WAIT1:$WAIT8 --job-name=$JOBNAME9a --time=$CALLVARIANT_TIME --cpus-per-task=$CALLVARIANT_CPU --mem=$CALLVARIANT_MEM \
            --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME9a".txt \
            --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,CALLER=$CALLER,i=$i \
            $IVDP_DIR/program/01.general_steps/slurm/9a.slurmCallingGVCFs.sb)
            echo $JOBID9a >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9a.txt
            #done
         done
         wait
         
      fi
   done
   
   WAIT9a=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9a.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')

fi
wait

            
if [[ $STEP_MDUP == "yes" && $STEP_BQSR == "no" && $STEP_VAR_CALL == "yes" ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9a.txt > /dev/null 2>&1
   
   for CALLER in $CALLER_PROGRAM; do
      ### VCF MODE
      if [[ $CALLER != "gatk4CombineGVCFs" && $CALLER != "gatk4GenomicsDBImport" && $CALLER != "gatk3" ]]; then
         for CHR in $(cat $CHROM_NAME); do
            JOBNAME9a="${ANALYSIS_ID}"_09.CallVCFs_"${CALLER}"_"${CHR}"
         
            JOBID9a=$(sbatch --dependency=afterok:$WAIT1:$WAIT7 --job-name=$JOBNAME9a --time=$CALLVARIANT_TIME --cpus-per-task=$CALLVARIANT_CPU --mem=$CALLVARIANT_MEM \
            --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME9a".txt \
            --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,CALLER=$CALLER,chr=$CHR \
            $IVDP_DIR/program/01.general_steps/slurm/9a.slurmCallingVCFs.sb)
            echo $JOBID9a >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9a.txt
         done
         wait
      ### CALL GVCF MODE
      else   
         for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
            i=$( echo $sample | awk -F\\ '{print $1}')
            
            #for CHR in $(cat $CHROM_NAME); do
            JOBNAME9a="${ANALYSIS_ID}"_09.CallGVCFs_"${CALLER}"_"${i}"
            
            JOBID9a=$(sbatch --dependency=afterok:$WAIT1:$WAIT7 --job-name=$JOBNAME9a --time=$CALLVARIANT_TIME --cpus-per-task=$CALLVARIANT_CPU --mem=$CALLVARIANT_MEM \
            --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME9a".txt \
            --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,CALLER=$CALLER,i=$i \
            $IVDP_DIR/program/01.general_steps/slurm/9a.slurmCallingGVCFs.sb)
            echo $JOBID9a >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9a.txt
            #done
         done
         wait
         
      fi
   done
   
   WAIT9a=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9a.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')

fi
wait


if [[ $STEP_MDUP == "no" && $STEP_BQSR == "no" && $STEP_VAR_CALL == "yes" ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9a.txt > /dev/null 2>&1
   
   for CALLER in $CALLER_PROGRAM; do
      ### VCF MODE
      if [[ $CALLER != "gatk4CombineGVCFs" && $CALLER != "gatk4GenomicsDBImport" && $CALLER != "gatk3" ]]; then
         for CHR in $(cat $CHROM_NAME); do
            JOBNAME9a="${ANALYSIS_ID}"_09.CallVCFs_"${CALLER}"_"${CHR}"
         
            JOBID9a=$(sbatch --dependency=afterok:$WAIT1 --job-name=$JOBNAME9a --time=$CALLVARIANT_TIME --cpus-per-task=$CALLVARIANT_CPU --mem=$CALLVARIANT_MEM \
            --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME9a".txt \
            --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,CALLER=$CALLER,chr=$CHR \
            $IVDP_DIR/program/01.general_steps/slurm/9a.slurmCallingVCFs.sb)
            echo $JOBID9a >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9a.txt
         done
         wait
      ### CALL GVCF MODE
      else   
         for sample in $(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2); do
            i=$( echo $sample | awk -F\\ '{print $1}')
            
            #for CHR in $(cat $CHROM_NAME); do
            JOBNAME9a="${ANALYSIS_ID}"_09.CallGVCFs_"${CALLER}"_"${i}"
            
            JOBID9a=$(sbatch --dependency=afterok:$WAIT1 --job-name=$JOBNAME9a --time=$CALLVARIANT_TIME --cpus-per-task=$CALLVARIANT_CPU --mem=$CALLVARIANT_MEM \
            --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME9a".txt \
            --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,CALLER=$CALLER,i=$i \
            $IVDP_DIR/program/01.general_steps/slurm/9a.slurmCallingGVCFs.sb)
            echo $JOBID9a >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9a.txt
            #done
         done
         wait
         
      fi
   done
   
   WAIT9a=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9a.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')

fi
wait


#----------------------------------#
##### STEP9b: CONCATENATE GVCF #####
#----------------------------------#
if [[ $STEP_VAR_CALL == "yes" ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9b.txt > /dev/null 2>&1

   for CALLER in $CALLER_PROGRAM; do
      if [[ $CALLER == "gatk4CombineGVCFs" || $CALLER == "gatk4GenomicsDBImport" || $CALLER == "gatk3" ]]; then
         for ALIGNER in $ALIGNER_PROGRAM; do
            JOBNAME9b="${ANALYSIS_ID}"_09.ConcatGVCF_"${ALIGNER}"_"${CALLER}"
            
            JOBID9b=$(sbatch --dependency=afterok:$WAIT1:$WAIT9a --job-name=$JOBNAME9b \
            --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME9b".txt \
            --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,ALIGNER=$ALIGNER,CALLER=$CALLER \
            $IVDP_DIR/program/01.general_steps/slurm/9b.slurmConcatGVCFs.sb)
            echo $JOBID9b >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9b.txt
         done
      fi
   done
   wait
   
   WAIT9b=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9b.txt 2>/dev/null | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
fi
wait


#-----------------------------#
##### STEP9c: GVCF TO VCF #####
#-----------------------------#
if [[ $STEP_VAR_CALL == "yes" && $STEP_GVCF_TO_VCF == "yes" ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9c.txt > /dev/null 2>&1

   for CALLER in $CALLER_PROGRAM; do   
      if [[ $CALLER == "gatk4CombineGVCFs" || $CALLER == "gatk4GenomicsDBImport" || $CALLER == "gatk3" ]]; then
         for ALIGNER in $ALIGNER_PROGRAM; do
            for CHR in $(cat $CHROM_NAME); do
               JOBNAME9c="${ANALYSIS_ID}"_09.GVCFtoVCF_"${ALIGNER}"_"${CALLER}"_"${CHR}"

               JOBID9c=$(sbatch --dependency=afterok:$WAIT1:$WAIT9b --job-name=$JOBNAME9c --time=$GVCF_TO_VCF_TIME --cpus-per-task=$GVCF_TO_VCF_CPU --mem=$GVCF_TO_VCF_MEM \
               --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME9c".txt \
               --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,ALIGNER=$ALIGNER,CALLER=$CALLER,chr=$CHR \
               $IVDP_DIR/program/01.general_steps/slurm/9c.slurmGVCF_TO_VCF.sb)
               echo $JOBID9c >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9c.txt
            done
         done
      fi
   done
   wait
      
   WAIT9c=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9c.txt 2>/dev/null | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')

fi
wait


if [[ $STEP_VAR_CALL == "no" && $STEP_GVCF_TO_VCF == "yes" ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9c.txt > /dev/null 2>&1

   for CALLER in $CALLER_PROGRAM; do   
      if [[ $CALLER == "gatk4CombineGVCFs" || $CALLER == "gatk4GenomicsDBImport" || $CALLER == "gatk3" ]]; then
         for ALIGNER in $ALIGNER_PROGRAM; do
            for CHR in $(cat $CHROM_NAME); do
               JOBNAME9c="${ANALYSIS_ID}"_09.GVCFtoVCF_"${ALIGNER}"_"${CALLER}"_"${CHR}"

               JOBID9c=$(sbatch --dependency=afterok:$WAIT1 --job-name=$JOBNAME9c --time=$GVCF_TO_VCF_TIME --cpus-per-task=$GVCF_TO_VCF_CPU --mem=$GVCF_TO_VCF_MEM \
               --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME9c".txt \
               --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,ALIGNER=$ALIGNER,CALLER=$CALLER,chr=$CHR \
               $IVDP_DIR/program/01.general_steps/slurm/9c.slurmGVCF_TO_VCF.sb)
               echo $JOBID9c >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9c.txt
            done
         done
      fi
   done
   wait
      
   WAIT9c=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9c.txt 2>/dev/null | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')

fi
wait


#-----------------------------------------#
##### STEP9d: CONCATENATE CHR OF VCFs #####
#-----------------------------------------#
if [[ $WAIT9a && $WAIT9c ]]; then
WAIT9_JOIN_VCF=$(echo $WAIT9a:$WAIT9c)
fi

if [[ $WAIT9a && -z $WAIT9c ]]; then
WAIT9_JOIN_VCF=$(echo $WAIT9a)
fi

if [[ -z $WAIT9a && $WAIT9c ]]; then
WAIT9_JOIN_VCF=$(echo $WAIT9c)
fi
wait


### Merge vcf for all the callers or only for gvcf callers
if [[ $STEP_VAR_CALL == "yes" && $STEP_GVCF_TO_VCF == "yes" ]] || [[ $STEP_VAR_CALL == "no" && $STEP_GVCF_TO_VCF == "yes" ]]; then
   rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9d.txt > /dev/null 2>&1
   JOBNAME9d="${ANALYSIS_ID}"_09.ConcatVCF

   JOBID9d=$(sbatch --dependency=afterok:$WAIT1:$WAIT9_JOIN_VCF --job-name=$JOBNAME9d \
   --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME9d".txt \
   --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME \
   $IVDP_DIR/program/01.general_steps/slurm/9d.slurmConcatVCFs.sb)
   echo $JOBID9d >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9d.txt
   
   WAIT9=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9d.txt 2>/dev/null  | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')

fi
wait


### Merge vcf only for non gvcf callers
# grep -vw: v means pattern not in strings; w means exact pattern as string; q means not to print grep result
if [[ $STEP_VAR_CALL == "yes" && $STEP_GVCF_TO_VCF == "no" ]]; then
   if echo $CALLER_PROGRAM | grep -wvq "gatk4CombineGVCFs\|gatk4GenomicsDBImport\|gatk3"; then
      rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9d.txt > /dev/null 2>&1
      JOBNAME9d="${ANALYSIS_ID}"_09.ConcatVCF
   
      JOBID9d=$(sbatch --dependency=afterok:$WAIT1:$WAIT9_JOIN_VCF --job-name=$JOBNAME9d \
      --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME9d".txt \
      --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME \
      $IVDP_DIR/program/01.general_steps/slurm/9d.slurmConcatVCFs.sb)
      echo $JOBID9d >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9d.txt
      
      WAIT9=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid9d.txt 2>/dev/null  | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
   fi
fi
wait


#----------------------------#
##### STEP10: FILTER VCF #####
#----------------------------#
if [[ $STEP_VCF_FILTER == "yes" ]]; then
   if [[ $STEP_VAR_CALL == "yes" || $STEP_GVCF_TO_VCF == "yes" ]]; then
      rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid10.txt > /dev/null 2>&1
      JOBNAME10="${ANALYSIS_ID}"_10.FilterVCFs

      JOBID10=$(sbatch --dependency=afterok:$WAIT1:$WAIT9 --job-name=$JOBNAME10 --time=$FILTERVCF_TIME --cpus-per-task=$FILTERVCF_CPU --mem=$FILTERVCF_MEM \
      --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME10".txt \
      --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME \
      $IVDP_DIR/program/01.general_steps/slurm/10.slurmFilterVCF.sb)
      echo $JOBID10 >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid10.txt
   
      WAIT10=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid10.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
   
   fi
fi
wait


if [[ $STEP_VCF_FILTER == "yes" ]]; then
   if [[ $STEP_VAR_CALL == "no" && $STEP_GVCF_TO_VCF == "no" ]]; then
      rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid10.txt > /dev/null 2>&1
      JOBNAME10="${ANALYSIS_ID}"_10.FilterVCFs

      JOBID10=$(sbatch --dependency=afterok:$WAIT1 --job-name=$JOBNAME10 --time=$FILTERVCF_TIME --cpus-per-task=$FILTERVCF_CPU --mem=$FILTERVCF_MEM \
      --output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"$JOBNAME10".txt \
      --export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME \
      $IVDP_DIR/program/01.general_steps/slurm/10.slurmFilterVCF.sb)
      echo $JOBID10 >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid10.txt
   
      WAIT10=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid10.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
   
   fi
fi
wait


#-----------------------------#
##### STEP11: HTML OUTPUT #####
#-----------------------------#
rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_wait_html.txt > /dev/null 2>&1
for job in $(ls $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_jobid*);do
   cat $job >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_wait_html.txt
done
wait
WAIT11=$(cat $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_wait_html.txt | awk '{print $4}' | tr '\n' ':' | sed 's/.$//')
rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_wait_html.txt > /dev/null 2>&1
JOBNAME11="${ANALYSIS_ID}"_"HTML"

sbatch --dependency=afterok:$WAIT11 --job-name=$JOBNAME11 \
--output=$OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/$JOBNAME11.txt \
--export=ANALYSIS_ID=$ANALYSIS_ID,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME \
$IVDP_DIR/program/01.general_steps/slurm/11.slurmHTML.sb


