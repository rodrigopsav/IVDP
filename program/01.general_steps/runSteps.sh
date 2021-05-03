#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate ivdp
conda activate --stack r-env4.0
export CONDA_DEFAULT_ENV=$CONDA_DEFAULT_ENV
wait


echo
echo
echo "#@#############################################################"
echo "#@                  NUMBER OF SAMPLES AND INDIVIDUALS "
echo "#@   Number of fastq files (samples): $N_SAMPLES"
echo "#@   Number of individuals: $N_INDIVIDUALS"
echo "#@"
echo "#@                         INPUT PARAMETERS "
echo "#@   INPUT FODLER: $INPUT_DIR"
echo "#@   OUTPUT FOLDER: $OUTPUT_DIR"
echo "#@   NAME OF ANALYSIS: $OUTPUT_NAME"
echo "#@   TYPE OF ANALYSIS: $(if [[ $ANALYSIS == 1 ]]; then echo "WGS"; else echo "RNA-Seq"; fi) "
echo "#@   TYPE OF FILES: $(if [[ $TYPE == "se" ]]; then echo "Single-end"; else echo "Paired-end"; fi) "
echo "#@   LIST OF SAMPLES: $LIST_SAMPLES"
echo "#@"
echo "#@                   REFERENCE GENOME PARAMETERS "
echo "#@   NAME OF REFERENCE GENOME: $ASSEMBLY"
echo "#@   REFERENCE GENOME FILE: $REFERENCE"
echo "#@   VARIANTS REFERENCE GENOME FILE: $VARIANTS"
echo "#@   ANNOTATION REFERENCE GENOME FILE: $ANNOTATION"
echo "#@"
echo "#@                         PIPELINE STEPS "
echo "#@   QUALITY CONTROL AND ADAPTER TRIMMING: $STEP_QC_TRIM"
echo "#@   ALIGMNMENT: $STEP_ALIGNMENT"
echo "#@   GENE COUNTS: $STEP_GENECOUNT"
echo "#@   MARK DUPLICATED READS: $STEP_MDUP"
echo "#@   BASE QUALITY SCORE RECALIBRATION: $STEP_BQSR"
echo "#@   VARIANT CALLING: $STEP_VAR_CALL"
echo "#@   CONVERT GVCF TO VCF: $STEP_GVCF_TO_VCF"
echo "#@   FILTER VCF: $STEP_VCF_FILTER"
echo "#@"
echo "#@                        OTHER PARAMETERS "
echo "#@   ALIGNER PROGRAMS: " $ALIGNER_PROGRAM
echo "#@   VARIANT CALLING PROGRAMS: " $CALLER_PROGRAM
echo "#@   GENE COUNT PROGRAM: " $GENECOUNT_PROGRAM
echo "#@   FEATURE TYPE FOR GENE COUNT: " $FEATURE_TYPE
echo "#@"
echo "#@   ADAPTERS FOR TRIMMING: $ADAPTER_TRIMMOMATIC"
echo "#@   MINIMUM READ LENGTH: $MIN_READ_LENGTH"
echo "#@   MAXIMUN READ LENGTH: $MAX_READ_LENGTH"
echo "#@   MINIMUN DEPTH FOR FILTERED LOCI: $MIN_DEPTH"
echo "#@   MINIMUN ALLELE FREQUENCY FOR FILTERED LOCI: $MAF"
echo "#@   MAXIMUN GENOTYPE MISSING RATE: $MISSING"
echo "#@   COMBINING VARIANTS FROM DIFFERENT VCFs: " $COMBINE_VCF
echo "#@"
if [[ "$DOWNSAMPLING_FASTQ" != 0 ]]; then echo "#@   DOWNSAMPLING COVERAGE (FASTQ): $DOWNSAMPLING_FASTQ"; fi
if [[ "$DOWNSAMPLING_BAM" != 0 ]]; then echo "#@   DOWNSAMPLING COVERAGE (BAM): $DOWNSAMPLING_BAM"; fi
echo "#@   NUMBER OF THREADS: $THREADS"
echo "#@   BATCH OF SAMPLES: $BATCH"
echo "#@   KEEP LOG FILES: $KEEP_LOG"
echo "#@   KEEP TRIMMED FASTQ AND BAM FILES: $KEEP_INTERMEDIATE_DATA" 
echo "#@"
D1=$(date "+%D    %T")
echo "#@    Date and Time: $D1"
echo "#@#############################################################"


START_IVDP=$(date +%s)

### INDEX REFERENCE
source $IVDP_DIR/program/01.general_steps/runGenomeIndex.sh
wait

### SUMMARY STATISTICS FASTQ FILES
source $IVDP_DIR/program/03.qc_trim/fastq_stats.sh
wait

### DOWNSAMPLING FASTQ FILES
if [[ $DOWNSAMPLING_FASTQ != 0 ]]; then
   source $IVDP_DIR/program/01.general_steps/runDownSamplingFastq.sh
fi
wait

### QUALITY CONTROL AND TRIM FASTQ FILES
if [[ $STEP_QC_TRIM == "yes" ]]; then
   source $IVDP_DIR/program/01.general_steps/runQcTrim.sh
fi
wait

### ALIGNMENT
if [[ $STEP_ALIGNMENT == "yes" ]]; then
   source $IVDP_DIR/program/01.general_steps/runAlignment.sh
fi
wait

### DOWNSAMPLING BAM FILES
if [[ $DOWNSAMPLING_BAM != 0 ]]; then
   source $IVDP_DIR/program/01.general_steps/runDownSamplingBam.sh
fi
wait

### GENE COUNT
if [[ $STEP_GENECOUNT == "yes" ]]; then
   source $IVDP_DIR/program/01.general_steps/runGeneCount.sh
fi
wait

### MARK DUPLICATED READS
if [[ $STEP_MDUP == "yes" ]]; then
   source $IVDP_DIR/program/01.general_steps/runAlignCorrectionMdupSplitncigar.sh
fi
wait

### BASE QUALITY SCORE RECALIBRATION
if [[ $STEP_BQSR == "yes" ]]; then
   if [[ "$BQSR_PROGRAM" == "bqsrgatk3" ]]; then
      conda activate --stack gatk3
      export CONDA_DEFAULT_ENV=$CONDA_DEFAULT_ENV
      source $IVDP_DIR/program/01.general_steps/runAlignCorrectionBQSR.sh
      wait
   else
      source $IVDP_DIR/program/01.general_steps/runAlignCorrectionBQSR.sh
      wait
   fi
fi
wait

### VARIANT CALLING 
if [[ $STEP_VAR_CALL == "yes" ]]; then
   for CALLER in $CALLER_PROGRAM; do
      export CALLER=$CALLER
      
      ### VCF MODE
      if [[ $CALLER != "gatk4CombineGVCFs" && $CALLER != "gatk4GenomicsDBImport" && $CALLER != "gatk3" && $CALLER != "freebayes" && $CALLER != "freebayes-parallel" && $CALLER != "platypus" ]]; then
         source $IVDP_DIR/program/01.general_steps/runCallingVCFs.sh
         wait
         
      elif [[ $CALLER == "freebayes" || $CALLER == "freebayes-parallel" || $CALLER == "platypus" ]]; then
         conda activate --stack ivdp2
         export CONDA_DEFAULT_ENV=$CONDA_DEFAULT_ENV
         source $IVDP_DIR/program/01.general_steps/runCallingVCFs.sh
         wait
         
      ### CALL GVCF MODE
      elif [[ $CALLER == "gatk3" ]]; then
         conda activate --stack gatk3
         export CONDA_DEFAULT_ENV=$CONDA_DEFAULT_ENV
         source $IVDP_DIR/program/01.general_steps/runCallingGVCFs.sh
         wait
         
      else  
         source $IVDP_DIR/program/01.general_steps/runCallingGVCFs.sh
         wait
      fi
      
   done
   wait

multiqc -f -n report_vcf_raw.html "$STATS_VCF_RAW"/*.stats -o "$REPORT_VCF_RAW" > /dev/null 2>&1
wait
fi
wait


### CONVERT GVCF TO VCF
if [[ $STEP_GVCF_TO_VCF == "yes" ]]; then
   for CALLER in $CALLER_PROGRAM; do
      export CALLER=$CALLER
      
      if [[ $CALLER == "gatk3" ]]; then
         conda activate --stack gatk3
         export CONDA_DEFAULT_ENV=$CONDA_DEFAULT_ENV     
         source $IVDP_DIR/program/01.general_steps/runCallingGVCF_TO_VCF.sh
      fi
      wait 
      
      
      if [[ $CALLER == "gatk4CombineGVCFs" || $CALLER == "gatk4GenomicsDBImport" ]]; then
         source $IVDP_DIR/program/01.general_steps/runCallingGVCF_TO_VCF.sh
      fi
      wait
      
   done
   wait

multiqc -f -n report_vcf_raw.html "$STATS_VCF_RAW"/*.stats -o "$REPORT_VCF_RAW" > /dev/null 2>&1
wait
fi
wait


### FILTER VCF
if [[ $STEP_VCF_FILTER == "yes" ]]; then
   source $IVDP_DIR/program/01.general_steps/runFilterVCF.sh
fi
wait


### CREATING HTML OUTPUT
if [[ $KEEP_INTERMEDIATE_DATA == "yes" ]]; then
   source $IVDP_DIR/program/01.general_steps/runOutputHtml1.sh > "$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/"$OUTPUT_NAME"_"$ASSEMBLY".html
else
   source $IVDP_DIR/program/01.general_steps/runOutputHtml2.sh > "$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/"$OUTPUT_NAME"_"$ASSEMBLY".html
fi
wait

source $IVDP_DIR/program/01.general_steps/runOutputHtml3.sh > "$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/report.html
cd "$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"
mkdir -p "$OUTPUT_NAME"_"$ASSEMBLY"_report 2> /dev/null
cp -r report "$OUTPUT_NAME"_"$ASSEMBLY"_report 2> /dev/null
cp -r report.html "$OUTPUT_NAME"_"$ASSEMBLY"_report 2> /dev/null
zip -q -u -r "$OUTPUT_NAME"_"$ASSEMBLY"_report.zip "$OUTPUT_NAME"_"$ASSEMBLY"_report 2> /dev/null
rm -r "$OUTPUT_NAME"_"$ASSEMBLY"_report report.html 2> /dev/null
wait
echo


END_IVDP=$(date +%s)
DIFF_IVDP=$(( $END_IVDP - $START_IVDP ))

echo
echo "@ WORKFLOW FLOW TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF_IVDP/3600)) $(($DIFF_IVDP%3600/60)) $(($DIFF_IVDP%60))) "
echo "@ Completed job "
echo "@#############################################################"
echo
echo