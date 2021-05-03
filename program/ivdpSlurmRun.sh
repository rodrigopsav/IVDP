#!/bin/bash

##########################################################
########### EXPORT VARIABLES OF PARAMETER FILE ###########
##########################################################

#
#
#

##########################################################
######################## INPUT DATA ######################

export INPUT_DIR=$(readlink -f $(cat $PARAMETERS | grep "|INPUT_DIR" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1) 2> /dev/null)
export OUTPUT_DIR=$(readlink -f $(cat $PARAMETERS | grep "|OUTPUT_DIR" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1) 2> /dev/null)
export OUTPUT_NAME=$(cat $PARAMETERS | grep "|OUTPUT_NAME" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export LIST_SAMPLES=$(cat $PARAMETERS | grep "|LIST_SAMPLES" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export ANALYSIS=$(cat $PARAMETERS | grep "|ANALYSIS" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export EXTENSION_SE=$(cat $PARAMETERS | grep "|EXTENSION_SE" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export EXTENSION_PE1=$(cat $PARAMETERS | grep "|EXTENSION_PE1" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export EXTENSION_PE2=$(cat $PARAMETERS | grep "|EXTENSION_PE2" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)


##########################################################
####################### INPUT GENOME #####################

export ASSEMBLY=$(cat $PARAMETERS | grep "|ASSEMBLY" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export REFERENCE=$(readlink -f $(cat $PARAMETERS | grep "|REFERENCE" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1) 2> /dev/null)
export VARIANTS=$(readlink -f $(cat $PARAMETERS | grep "|VARIANTS" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1) 2> /dev/null)
export ANNOTATION=$(readlink -f $(cat $PARAMETERS | grep "|ANNOTATION" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1) 2> /dev/null)


##########################################################
####################### WORKFLOW STEPS ###################

export STEP_QC_TRIM=$(cat $PARAMETERS | grep "|STEP_QC_TRIM" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export STEP_ALIGNMENT=$(cat $PARAMETERS | grep "|STEP_ALIGNMENT" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export STEP_GENECOUNT=$(cat $PARAMETERS | grep "|STEP_GENECOUNT" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export STEP_MDUP=$(cat $PARAMETERS | grep "|STEP_MDUP" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export STEP_BQSR=$(cat $PARAMETERS | grep "|STEP_BQSR" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export STEP_VAR_CALL=$(cat $PARAMETERS | grep "|STEP_VAR_CALL" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export STEP_GVCF_TO_VCF=$(cat $PARAMETERS | grep "|STEP_GVCF_TO_VCF" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export STEP_VCF_FILTER=$(cat $PARAMETERS | grep "|STEP_VCF_FILTER" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)


##########################################################
#################### CHOOSE THE ANALYSIS #################
#################### CHOOSE THE PROGRAMS #################

export ALIGNER_PROGRAM=$(cat $PARAMETERS | grep "|ALIGNER_PROGRAM" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1 | tr "," "\n")
export CALLER_PROGRAM=$(cat $PARAMETERS | grep "|CALLER_PROGRAM" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1 | tr "," "\n")
export GENECOUNT_PROGRAM=$(cat $PARAMETERS | grep "|GENECOUNT_PROGRAM" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1 | tr "," "\n")
export FEATURE_TYPE=$(cat $PARAMETERS | grep "|FEATURE_TYPE" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1 | tr "," "\n")
export MDUP_PROGRAM=$(cat $PARAMETERS | grep "|MDUP_PROGRAM" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1 | tr "," "\n")
export BQSR_PROGRAM=$(cat $PARAMETERS | grep "|BQSR_PROGRAM" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1 | tr "," "\n")


##########################################################
#################### CHOOSE THE ANALYSIS #################
##################### GENERAL PARAMETERS #################

export STATS_ALIGNMENT=$(cat $PARAMETERS | grep "|STATS_ALIGNMENT" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
#export CALL_BY_CHROM=$(cat $PARAMETERS | grep "|CALL_BY_CHROM" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export BP_BY_CHROM=$(cat $PARAMETERS | grep "|BP_BY_CHROM" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export CHROM_SET=$(cat $PARAMETERS | grep "|CHROM_SET" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export SELECT_CHROM=$(cat $PARAMETERS | grep "|SELECT_CHROM" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1 | tr "," "\n")
export INCLUDE_CHROM_X=$(cat $PARAMETERS | grep "|INCLUDE_CHROM_X" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export INCLUDE_CHROM_Y=$(cat $PARAMETERS | grep "|INCLUDE_CHROM_Y" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export INCLUDE_CHROM_MT=$(cat $PARAMETERS | grep "|INCLUDE_CHROM_MT" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export ADAPTER_TRIMMOMATIC=$(cat $PARAMETERS | grep "|ADAPTER_TRIMMOMATIC" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export MIN_READ_LENGTH=$(cat $PARAMETERS | grep "|MIN_READ_LENGTH" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export MAX_READ_LENGTH=$(cat $PARAMETERS | grep "|MAX_READ_LENGTH" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export DOWNSAMPLING_FASTQ=$(cat $PARAMETERS | grep "|DOWNSAMPLING_FASTQ" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export DOWNSAMPLING_BAM=$(cat $PARAMETERS | grep "|DOWNSAMPLING_BAM" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export MIN_DEPTH=$(cat $PARAMETERS | grep "|MIN_DEPTH" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export MAF=$(cat $PARAMETERS | grep "|MAF" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export MISSING=$(cat $PARAMETERS | grep "|MISSING" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export COMBINE_VCF=$(cat $PARAMETERS | grep "|COMBINE_VCF" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1 | tr "," "\n")
#export THREADS=$(cat $PARAMETERS | grep "|THREADS" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
#export BATCH=$(cat $PARAMETERS | grep "|BATCH" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export KEEP_LOG=$(cat $PARAMETERS | grep "|KEEP_LOG" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export KEEP_INTERMEDIATE_DATA=$(cat $PARAMETERS | grep "|KEEP_INTERMEDIATE_DATA" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)

##########################################################

#
#
#

### CHECK CONDA ENVIRONMENTS
eval "$(conda shell.bash hook)"
conda activate ivdp
conda activate --stack r-env4.0
if [[ "$CONDA_DEFAULT_ENV" != "ivdp" && "$CONDA_DEFAULT_ENV" != "r-env4.0" ]]; then
   echo
   echo "ERROR: conda ivdp environment was not found. Check if the IVDP dependencies are installed (install_ivdp_dependencies/install_ivdp_dependencies.sh)"
   echo
   exit 1
fi
wait


conda activate --stack ivdp2
if [[ "$CONDA_DEFAULT_ENV" != "ivdp2" ]]; then
   echo
   echo "ERROR: conda ivdp2 environment was not found. Check if the IVDP dependencies are installed (install_ivdp_dependencies/install_ivdp_dependencies.sh)"
   echo
   exit 1
fi
wait


source $IVDP_DIR/program/01.general_steps/runCheckConda.sh
wait
conda deactivate
wait
#export CONDA_DEFAULT_ENV=$CONDA_DEFAULT_ENV
wait


if $(echo $CALLER_PROGRAM | grep -q "gatk3"); then
   conda activate --stack gatk3
   if [[ "$CONDA_DEFAULT_ENV" != "gatk3" ]]; then
      echo
      echo "ERROR: conda gatk3 environment was not found. Check if gatk3 module is installed"
      echo
      exit 1
   fi
   conda deactivate
   wait
fi
wait

### CHECK IVDP PARAMETERS
source $IVDP_DIR/program/01.general_steps/runCheckParameters.sh
wait

### CREATE LIST OF SAMPLES
export ANALYSIS_ID=$(echo $RANDOM)
source $IVDP_DIR/program/01.general_steps/runCreateListSamples.sh
wait
sed -i 's/\t/\\/g' ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt
dos2unix ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt 2> /dev/null
wait

### CREATE CHROMOSOME NAMES
source $IVDP_DIR/program/01.general_steps/runChromNames.sh
wait

##########################################################
############## CREATE AUXILIARY PARAMETER FILE ###########
##########################################################
mkdir -p $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm > /dev/null 2>&1
rm $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt > /dev/null 2>&1


echo "|IVDP_DIR="$IVDP_DIR > $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|INPUT_DIR="$INPUT_DIR >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|OUTPUT_DIR="$OUTPUT_DIR >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|OUTPUT_NAME="$OUTPUT_NAME >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|LIST_SAMPLES="$LIST_SAMPLES >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|ANALYSIS="$ANALYSIS >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|EXTENSION_SE="$EXTENSION_SE >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|EXTENSION_PE1="$EXTENSION_PE1 >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|EXTENSION_PE2="$EXTENSION_PE2 >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|ASSEMBLY="$ASSEMBLY >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|REFERENCE="$REFERENCE >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|VARIANTS="$VARIANTS >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|ANNOTATION="$ANNOTATION >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|STEP_QC_TRIM="$STEP_QC_TRIM >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|STEP_ALIGNMENT="$STEP_ALIGNMENT >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|STEP_GENECOUNT="$STEP_GENECOUNT >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|STEP_MDUP="$STEP_MDUP >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|STEP_BQSR="$STEP_BQSR >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|STEP_VAR_CALL="$STEP_VAR_CALL >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|STEP_GVCF_TO_VCF="$STEP_GVCF_TO_VCF >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|STEP_VCF_FILTER="$STEP_VCF_FILTER >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|MDUP_PROGRAM="$MDUP_PROGRAM >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|BQSR_PROGRAM="$BQSR_PROGRAM >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|STATS_ALIGNMENT="$STATS_ALIGNMENT >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|CALL_BY_CHROM="$CALL_BY_CHROM >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|BP_BY_CHROM="$BP_BY_CHROM >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|CHROM_SET="$CHROM_SET >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|INCLUDE_CHROM_X="$INCLUDE_CHROM_X >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|INCLUDE_CHROM_Y="$INCLUDE_CHROM_Y >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|INCLUDE_CHROM_MT="$INCLUDE_CHROM_MT >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|ADAPTER_TRIMMOMATIC="$ADAPTER_TRIMMOMATIC >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|MIN_READ_LENGTH="$MIN_READ_LENGTH >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|MAX_READ_LENGTH="$MAX_READ_LENGTH >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|DOWNSAMPLING_FASTQ="$DOWNSAMPLING_FASTQ >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|DOWNSAMPLING_BAM="$DOWNSAMPLING_BAM >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|MIN_DEPTH="$MIN_DEPTH >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|MAF="$MAF >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|MISSING="$MISSING >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|BATCH="$BATCH >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|KEEP_LOG="$KEEP_LOG >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|KEEP_INTERMEDIATE_DATA="$KEEP_INTERMEDIATE_DATA >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|TYPE="$TYPE >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|PARAMETERS="$PARAMETERS >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo
echo "|CHROM_NAME="$CHROM_NAME >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|CHROM_SIZE="$CHROM_SIZE >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|CHROM_POSITIONS="$CHROM_POSITIONS >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|ALIGNER_PROGRAM="$ALIGNER_PROGRAM | sed -e 's/\s\+/,/g' >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|CALLER_PROGRAM="$CALLER_PROGRAM | sed -e 's/\s\+/,/g' >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|GENECOUNT_PROGRAM="$GENECOUNT_PROGRAM | sed -e 's/\s\+/,/g' >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|FEATURE_TYPE="$FEATURE_TYPE | sed -e 's/\s\+/,/g' >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt
echo "|COMBINE_VCF="$COMBINE_VCF | sed -e 's/\s\+/,/g' >> $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm/"${ANALYSIS_ID}"_parAuxiliar.txt

### DEACTIVATE CONDA ENVIRONMENT
# This is very important, otherwise it'll conflict with conda env of runStepsSlurm.sh
# Remember that it is not possible to activate conda in subshells. So, all conda envs will be activated in from 1.sb to 13.sb files.
conda deactivate
conda deactivate

### RUN IVDP
mkdir -p $OUTPUT_DIR/ivdp_"$OUTPUT_NAME"/logSlurm > /dev/null 2>&1
source $IVDP_DIR/program/01.general_steps/slurm/runSlurmSteps.sh
echo "Running IVDP for HPCC with Slurm"
echo "ANALYSIS NAME: $OUTPUT_NAME"
echo "ANALYSIS ID: $ANALYSIS_ID"
echo "To check all the slurm jobs, run: $IVDP_DIR/ivdpLog.sh -l 2 -a $ANALYSIS_ID -o $OUTPUT_DIR/ivdp_${OUTPUT_NAME}"
echo "To kill this IVDP analysis, run:"
echo "    for job in \$(squeue -u $USER | grep "$ANALYSIS_ID" | awk '{print \$1}'); do scancel \$job; done"

