#!/bin/bash

echo
echo
echo "#@############################################################"
echo "#@                      CHECK IVDP PARAMETERS"
#echo


#############################################
########## LOAD DEFAULT PARAMETERS ##########
#############################################
COUNT_ERROR=0
COUNT_WARNING=0

################################
##### INPUT DATA VARIABLES #####
################################

### INPUT_DIR
if [[ -z "$INPUT_DIR" ]]; then
   echo "#@ ERROR: Missing input file directory in parameter file (INPUT_DIR variable) "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
else
   if [[ ! -d "$INPUT_DIR" ]]; then
      echo "#@ ERROR: Invalid input file directory in parameter file (INPUT_DIR variable) "
      COUNT_ERROR=$((COUNT_ERROR + 1))
      echo
   fi
fi

### OUTPUT_DIR
if [[ -z "$OUTPUT_DIR" ]]; then
   echo "#@ WARNING: Missing output directory in parameter file (OUTPUT_DIR variable) "
   echo "#@ The output folder will be created in $"HOME" directory"
   export OUTPUT_DIR=$HOME
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
fi

### OUTPUT_NAME
if [[ -z "$OUTPUT_NAME" ]]; then
   echo "#@ WARNING: Missing output name of the analysis in parameter file (OUTPUT_NAME variable) "
   echo "#@ The output name will be $"USER" "
   export OUTPUT_NAME=$USER
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
fi

### LIST_SAMPLES
if [[ -z "$LIST_SAMPLES" ]]; then
   export LIST_SAMPLES=none
elif [[ "${LIST_SAMPLES,,}" == "none" ]]; then
   export LIST_SAMPLES=none
else
   export LIST_SAMPLES=$(readlink -f $LIST_SAMPLES)
   if [[ ! -f "$LIST_SAMPLES" ]]; then
      echo "#@ ERROR: LIST_SAMPLES files does not exist. Please check its name in parameter file (LIST_SAMPLES variable) "
      COUNT_ERROR=$((COUNT_ERROR + 1))
      echo      
   fi
fi

### ANALYSIS
if [[ -z "$ANALYSIS" ]]; then
   echo "#@ ERROR: Missing analysis type in parameter file (ANALYSIS variable) "
   echo "#@ Choose ANALYSIS=1 (WGS) or ANALYSIS=2 (RNAseq)"
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

if [[ $ANALYSIS != "1" && $ANALYSIS != "2" ]]; then
   echo "#@ ERROR: Invalid option for analysis type in parameter file (ANALYSIS variable) "
   echo "#@ Choose ANALYSIS=1 (WGS) or ANALYSIS=2 (RNAseq)"
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

### EXTENSION and TYPE
if [[ -z "$EXTENSION_SE" && -z "$EXTENSION_PE1" && -z "$EXTENSION_PE2" ]]; then
   echo "#@ ERROR: Missing extension of fastq files in parameter file "
   echo "#@ Use: EXTENSION_SE or EXTENSION_PE1 and EXTENSION_PE2 variables"
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

if [[ "$EXTENSION_SE" && "$EXTENSION_PE1" && "$EXTENSION_PE2"  ]]; then
   echo "#@ ERROR: Impossible to detect if fastq files are single-end or paired-end "
   echo "#@ Use: EXTENSION_SE or EXTENSION_PE1 and EXTENSION_PE2 variables"
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

if [[ "$EXTENSION_SE" && "$EXTENSION_PE1" && -z "$EXTENSION_PE2" ]]; then
   echo "#@ ERROR: Impossible to detect if fastq files are single-end or paired-end "
   echo "#@ Use: EXTENSION_SE or EXTENSION_PE1 and EXTENSION_PE2 variables"
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

if [[ "$EXTENSION_SE" && -z "$EXTENSION_PE1" && "$EXTENSION_PE2" ]]; then
   echo "#@ ERROR: Impossible to detect if fastq files are single-end or paired-end "
   echo "#@ Use: EXTENSION_SE or EXTENSION_PE1 and EXTENSION_PE2 variables"
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

if [[ -z "$EXTENSION_SE" && "$EXTENSION_PE1" && -z "$EXTENSION_PE2" ]]; then
   echo "#@ ERROR: Specify both EXTENSION_PE1 and EXTENSION_PE2 "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

if [[ -z "$EXTENSION_SE" && -z "$EXTENSION_PE1" && "$EXTENSION_PE2" ]]; then
   echo "#@ ERROR: Specify both EXTENSION_PE1 and EXTENSION_PE2 "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

if [[ "$EXTENSION_SE" && -z "$EXTENSION_PE1" && -z "$EXTENSION_PE2" ]]; then
   export TYPE=se
fi

if [[ -z "$EXTENSION_SE" && "$EXTENSION_PE1" && "$EXTENSION_PE2" ]]; then
   export TYPE=pe
fi


##################################
##### INPUT GENOME VARIABLES #####
##################################

### ASSEMBLY
if [[ -z "$ASSEMBLY" ]]; then
   echo "#@ WARNING: Missing assembly name in parameter file (ASSEMBLY variable) "
   echo "#@ The ASSEMBLY name will be $"USER"ASSEMBLY "
   export ASSEMBLY="$USER"ASSEMBLY
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
fi

### REFERENCE
if [[ -z "$REFERENCE" ]]; then
   echo "#@ ERROR: Missing Reference.fa file directory in parameter file (REFERENCE variable) "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
else
   if [[ ! -d "${REFERENCE%/*}" ]]; then
      echo "#@ ERROR: Invalid Reference.fa directory in parameter file (REFERENCE variable) "
      COUNT_ERROR=$((COUNT_ERROR + 1))
      echo
   fi
fi

### VARIANTS
if [[ -z "$VARIANTS" || "${VARIANTS,,}" == "none" ]]; then
   echo "#@ WARNING: Missing Variants.vcf.gz file directory in parameter file (VARIANTS variable) "
   echo "#@ Base quality score recalibration step will be skipped"
   export VARIANTS=none
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
else
   if [[ ! -d "${VARIANTS%/*}" ]]; then
      echo "#@ ERROR: Invalid Variants.vcf.gz directory in parameter file (VARIANTS variable) "
      COUNT_ERROR=$((COUNT_ERROR + 1))
      echo
   fi
fi

### ANNOTATION
if [[ -z "$ANNOTATION" || "${ANNOTATION,,}" == "none" ]]; then
   echo "#@ WARNING: Missing Annotation.gtf file directory in parameter file (ANNOTATION variable) "
   echo "#@ Gene count step will be skipped in RNAseq analysis"
   echo "#@ VCF files will not be annotated"
   export ANNOTATION=none
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
else   
   if [[ ! -d "${ANNOTATION%/*}" ]]; then
      echo "#@ ERROR: Invalid Annotation.gtf directory in parameter file (ANNOTATION variable) "
      COUNT_ERROR=$((COUNT_ERROR + 1))
      echo
   fi
fi


###################################
##### WORKFLOW STEP VARIABLES #####
###################################

### STEP_QC_TRIM
# Defaults
if [[ -z "$STEP_QC_TRIM" ]]; then 
   export STEP_QC_TRIM=yes
fi

if [[ "${STEP_QC_TRIM,,}" != "no" && "${STEP_QC_TRIM,,}" != "yes" ]]; then
   echo "#@ WARNING: Invalid option in parameter file (STEP_QC_TRIM variable) "
   echo "#@ IVDP will consider STEP_QC_TRIM=yes"
   export STEP_QC_TRIM=yes
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
else
   export STEP_QC_TRIM=${STEP_QC_TRIM,,}
fi


### STEP_ALIGNMENT
# Defaults
if [[ -z "$STEP_ALIGNMENT" ]]; then 
   export STEP_ALIGNMENT=yes
fi

if [[ "${STEP_ALIGNMENT,,}" != "no" && "${STEP_ALIGNMENT,,}" != "yes" ]]; then
   echo "#@ WARNING: Invalid option in parameter file (STEP_ALIGNMENT variable) "
   echo "#@ IVDP will consider STEP_ALIGNMENT=yes"
   export STEP_ALIGNMENT=yes
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
else
   export STEP_ALIGNMENT=${STEP_ALIGNMENT,,}
fi


### STEP_GENECOUNT
# Defaults

if [[ "${ANNOTATION,,}" == "none" ]]; then
   export STEP_GENECOUNT=no
else
   if [[ -z "$STEP_GENECOUNT" && "$ANALYSIS" == "1" ]]; then 
      export STEP_GENECOUNT=no

   elif [[ "$STEP_GENECOUNT" && "$ANALYSIS" == "1" ]]; then 
      export STEP_GENECOUNT=no

   elif [[ -z "$STEP_GENECOUNT" && "$ANALYSIS" == "2" ]]; then 
      export STEP_GENECOUNT=yes

   else
   
      if [[ "${STEP_GENECOUNT,,}" != "no" && "${STEP_GENECOUNT,,}" != "yes" ]]; then
         echo "#@ WARNING: Invalid option in parameter file (STEP_GENECOUNT variable) "
         echo "#@ IVDP will consider STEP_GENECOUNT=yes"
         export STEP_GENECOUNT=yes
         COUNT_WARNING=$((COUNT_WARNING + 1))
         echo
      else
         export STEP_GENECOUNT=${STEP_GENECOUNT,,}
      fi
   fi
fi


### STEP_MDUP
# Defaults
if [[ -z "$STEP_MDUP" ]]; then 
   export STEP_MDUP=yes
fi

if [[ "${STEP_MDUP,,}" != "no" && "${STEP_MDUP,,}" != "yes" ]]; then
   echo "#@ WARNING: Invalid option in parameter file (STEP_MDUP variable) "
   echo "#@ IVDP will consider STEP_MDUP=yes"
   export STEP_MDUP=yes
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
else
   export STEP_MDUP=${STEP_MDUP,,}
fi


### STEP_BQSR
# Defaults
if [[ "${VARIANTS,,}" == "none" ]]; then
   export STEP_BQSR=no
else
   if [[ -z "$STEP_BQSR" ]]; then 
      export STEP_BQSR=yes
   fi
fi

if [[ "${STEP_BQSR,,}" != "no" && "${STEP_BQSR,,}" != "yes" ]]; then
   echo "#@ WARNING: Invalid option in parameter file variable (STEP_BQSR variable) "
   echo "#@ IVDP will consider STEP_BQSR=yes"
   export STEP_BQSR=yes
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
else
   export STEP_BQSR=${STEP_BQSR,,}
fi


### STEP_VAR_CALL
# Defaults
if [[ -z "$STEP_VAR_CALL" ]]; then 
   export STEP_VAR_CALL=yes
fi

if [[ "${STEP_VAR_CALL,,}" != "no" && "${STEP_VAR_CALL,,}" != "yes" ]]; then
   echo "#@ WARNING: Invalid option in parameter file (STEP_VAR_CALL variable) "
   echo "#@ IVDP will consider STEP_VAR_CALL=yes"
   export STEP_VAR_CALL=yes
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
else
   export STEP_VAR_CALL=${STEP_VAR_CALL,,}
fi


### STEP_GVCF_TO_VCF
# Defaults
if [[ -z "$STEP_GVCF_TO_VCF" ]]; then 
   export STEP_GVCF_TO_VCF=yes
fi

if [[ "${STEP_GVCF_TO_VCF,,}" != "no" && "${STEP_GVCF_TO_VCF,,}" != "yes" ]]; then
   echo "#@ WARNING: Invalid option in parameter file (STEP_GVCF_TO_VCF variable) "
   echo "#@ IVDP will consider STEP_GVCF_TO_VCF=yes"
   export STEP_GVCF_TO_VCF=yes
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
else
   export STEP_GVCF_TO_VCF=${STEP_GVCF_TO_VCF,,}
fi


### STEP_VCF_FILTER
# Defaults
if [[ -z "$STEP_VCF_FILTER" ]]; then 
   export STEP_VCF_FILTER=yes
fi

if [[ "${STEP_VCF_FILTER,,}" != "no" && "${STEP_VCF_FILTER,,}" != "yes" ]]; then
   echo "#@ WARNING: Invalid option in parameter file (STEP_VCF_FILTER variable) "
   echo "#@ IVDP will consider STEP_VCF_FILTER=yes"
   export STEP_VCF_FILTER=yes
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
else
   export STEP_VCF_FILTER=${STEP_VCF_FILTER,,}
fi


###################################
##### INPUT PROGRAM VARIABLES #####
###################################

### ALIGNER PROGRAMS
# Defaults
if [[ "$ANALYSIS" == 1 ]] && [[ -z "$ALIGNER_PROGRAM" ]]; then
   export ALIGNER_PROGRAM=bwa2
elif [[ "$ANALYSIS" == 2 ]] && [[ -z "$ALIGNER_PROGRAM" ]]; then
   export ALIGNER_PROGRAM=hisat2
else
   echo
fi

# Check
ALIGNER_SET="bbmap bowtie bowtie2 bwa bwa2 gsnap hisat2 star"
ALIGNER_PROGRAM=${ALIGNER_PROGRAM,,}

COUNT_PROGRAMS=0
COUNT_MATCH=0
for i in $ALIGNER_PROGRAM; do
   COUNT_PROGRAMS=$((COUNT_PROGRAMS + 1))
   for j in $ALIGNER_SET; do
   
      if [[ $i == $j ]]; then
         COUNT_MATCH=$((COUNT_MATCH + 1))
      fi
      
   done
done

if [[ "$COUNT_PROGRAMS" != "$COUNT_MATCH" ]]; then
   echo "#@ ERROR: Invalid aligner program in parameter file (ALIGNER_PROGRAM variable) "
   COUNT_ERROR=$((COUNT_ERROR + 1))
else
   export ALIGNER_PROGRAM=${ALIGNER_PROGRAM,,}
fi


### VARIANT CALLER PROGRAMS
# Defaults
if [[ "$ANALYSIS" == 1 ]] && [[ -z "$CALLER_PROGRAM" ]]; then
   export CALLER_PROGRAM=gatk4GenomicsDBImport
elif [[ "$ANALYSIS" == 2 ]] && [[ -z "$CALLER_PROGRAM" ]]; then
   export CALLER_PROGRAM=gatk4
else
   echo
fi

# Check
CALLER_SET="bcftools freebayes freebayes-parallel gatk4 gatk4CombineGVCFs gatk4GenomicsDBImport platypus varscan gatk3"
CALLER_PROGRAM=$CALLER_PROGRAM

COUNT_PROGRAMS=0
COUNT_MATCH=0
for i in $CALLER_PROGRAM; do
   COUNT_PROGRAMS=$((COUNT_PROGRAMS + 1))
   for j in $CALLER_SET; do

      if [[ $i == $j ]]; then
         COUNT_MATCH=$((COUNT_MATCH + 1))
      fi

   done
done

if [[ "$COUNT_PROGRAMS" != "$COUNT_MATCH" ]]; then
   echo "#@ ERROR: Invalid variant caller program in parameter file (CALLER_PROGRAM variable) "
   COUNT_ERROR=$((COUNT_ERROR + 1))
else
   export CALLER_PROGRAM=$CALLER_PROGRAM
fi


### GENE COUNT PROGRAMS
# Defaults
if [[ "$ANALYSIS" == 1 && "$STEP_GENECOUNT" == "no" ]]; then
   export GENECOUNT_PROGRAM=none
   export FEATURE_TYPE=none
fi

if [[ "$ANALYSIS" == 1 && "$STEP_GENECOUNT" == "yes" ]]; then
   export GENECOUNT_PROGRAM=none
   export FEATURE_TYPE=none
fi

if [[ "$ANALYSIS" == 2 && "$STEP_GENECOUNT" == "no" ]]; then
   export GENECOUNT_PROGRAM=none
   export FEATURE_TYPE=none
fi

if [[ "$ANALYSIS" == 2 && "$STEP_GENECOUNT" == "yes" ]]; then
   if [[ "$GENECOUNT_PROGRAM" ]]; then
   
      # Check
      GENECOUNT_SET="none htseq featurecounts"
      GENECOUNT_PROGRAM=${GENECOUNT_PROGRAM,,}
      
      COUNT_PROGRAMS=0
      COUNT_MATCH=0
      for i in $GENECOUNT_PROGRAM; do
         COUNT_PROGRAMS=$((COUNT_PROGRAMS + 1))
         for j in $GENECOUNT_SET; do
      
            if [[ $i == $j ]]; then
               COUNT_MATCH=$((COUNT_MATCH + 1))
            fi
      
         done
      done 

      if [[ "$COUNT_PROGRAMS" != "$COUNT_MATCH" ]]; then
         echo "#@ ERROR: Invalid gene count programs in parameter file (GENECOUNT_PROGRAM variable) "
         echo "#@ Choose GENECOUNT_PROGRAM=htseq or GENECOUNT_PROGRAM=featurecounts"
         COUNT_ERROR=$((COUNT_ERROR + 1))
      else
         export GENECOUNT_PROGRAM=${GENECOUNT_PROGRAM,,}
      fi
   
   elif [[ -z "$GENECOUNT_PROGRAM" ]]; then
      export GENECOUNT_PROGRAM=htseq
      export FEATURE_TYPE=exon
      
   else
      echo ""

   fi
fi


### FEATURE_TYPE
if [[ "$ANALYSIS" == 2 && "$STEP_GENECOUNT" == "yes" ]]; then
   if [[ -z "$FEATURE_TYPE" ]]; then
      export FEATURE_TYPE=exon
   else
      export FEATURE_TYPE=$FEATURE_TYPE
   fi
fi


### MARK DUPLICATED READS PROGRAM
# Defaults
if [[ -z "$MDUP_PROGRAM" ]]; then 
   export MDUP_PROGRAM=sambamba
fi

# Check
MDUP_SET="markdupspark picard sambamba"
MDUP_PROGRAM=${MDUP_PROGRAM,,}

COUNT_PROGRAMS=0
COUNT_MATCH=0
for i in $MDUP_PROGRAM; do
   COUNT_PROGRAMS=$((COUNT_PROGRAMS + 1))
   for j in $MDUP_SET; do

      if [[ $i == $j ]]; then
         COUNT_MATCH=$((COUNT_MATCH + 1))
      fi

   done
done

if [[ "$COUNT_PROGRAMS" != "$COUNT_MATCH" ]]; then
   echo "#@ ERROR: Invalid mark duplicated read programs in parameter file (MDUP_PROGRAM variable) "
   echo "#@ Choose MDUP_PROGRAM=sambamba or MDUP_PROGRAM=markdupspark or MDUP_PROGRAM=picard"
   COUNT_ERROR=$((COUNT_ERROR + 1))
else
   export MDUP_PROGRAM=${MDUP_PROGRAM,,}
fi


### BQSR PROGRAM
# Defaults
if [[ -z "$BQSR_PROGRAM" ]]; then 
   export BQSR_PROGRAM=bqsrspark
fi

# Check
BQSR_SET="bqsr bqsrspark bqsrgatk3"
BQSR_PROGRAM=${BQSR_PROGRAM,,}

COUNT_PROGRAMS=0
COUNT_MATCH=0
for i in $BQSR_PROGRAM; do
   COUNT_PROGRAMS=$((COUNT_PROGRAMS + 1))
   for j in $BQSR_SET; do

      if [[ $i == $j ]]; then
         COUNT_MATCH=$((COUNT_MATCH + 1))
      fi

   done
done

if [[ "$COUNT_PROGRAMS" != "$COUNT_MATCH" ]]; then
   echo "#@ ERROR: Invalid BQSR programs in parameter file (BQSR_PROGRAM variable) "
   echo "#@ Choose BQSR_PROGRAM=bqsrspark or BQSR_PROGRAM=bqsr"
   COUNT_ERROR=$((COUNT_ERROR + 1))
else
   export BQSR_PROGRAM=${BQSR_PROGRAM,,}
fi


##############################
##### GENERAL PARAMETERS #####
##############################

### STATS_ALIGNMENT
if [[ -z "$STATS_ALIGNMENT" ]]; then
   export STATS_ALIGNMENT=yes
fi

if [[ "${STATS_ALIGNMENT,,}" != "no" && "${STATS_ALIGNMENT,,}" != "yes" ]]; then
   echo "#@ WARNING: Invalid option for STATS_ALIGNMENT variable in parameter file "
   echo "#@ IVDP will consider STATS_ALIGNMENT=yes"
   export STATS_ALIGNMENT=yes
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
else
   export STATS_ALIGNMENT=${STATS_ALIGNMENT,,}
fi


### CALL_BY_CHROM
if [[ -z "$CALL_BY_CHROM" ]]; then
   export CALL_BY_CHROM=yes
fi

if [[ "${CALL_BY_CHROM,,}" != "no" && "${CALL_BY_CHROM,,}" != "yes" ]]; then
   echo "#@ WARNING: Invalid option for CALL_BY_CHROM variable in parameter file "
   echo "#@ IVDP will consider CALL_BY_CHROM=yes"
   export CALL_BY_CHROM=yes
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
else
   export CALL_BY_CHROM=${CALL_BY_CHROM,,}
fi


### BP_BY_CHROM
if [[ -z "$BP_BY_CHROM" ]]; then
   export BP_BY_CHROM=all
fi

if [[ "${BP_BY_CHROM,,}" != "all" ]]; then
   if [[ "${BP_BY_CHROM}" =~ ^[0-9]+$ ]]; then
      export BP_BY_CHROM=$BP_BY_CHROM
   else
      echo "#@ WARNING: Invalid option for BP_BY_CHROM variable in parameter file "
      echo "#@ IVDP will consider BP_BY_CHROM=all"
      export BP_BY_CHROM=all
      COUNT_WARNING=$((COUNT_WARNING + 1))
      echo
   fi
fi


### CHROM_SET
if [[ -z "$CHROM_SET" ]]; then
   echo "#@ ERROR: Missing value in parameter file (CHROM_SET variable) "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi


### SELECT_CHROM
# Defaults
if [[ -z "$SELECT_CHROM" ]]; then 
   export SELECT_CHROM=$(seq 1 1 $CHROM_SET)
else
   if [[ "$SELECT_CHROM" == "all" ]]; then
      export SELECT_CHROM=$(seq 1 1 $CHROM_SET)
   else
      export SELECT_CHROM=$SELECT_CHROM
   fi
fi


if [[ "$(echo $SELECT_CHROM | tr " " "\n" | wc -l)" -gt "$CHROM_SET" ]]; then
   echo "#@ ERROR: Number of selected autosomes (SELECT_CHROM variable) is greater than CHROM_SET."
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi


for chr in $SELECT_CHROM; do 
   if [[ "$chr" -gt "$CHROM_SET" ]]; then 
      echo "#@ ERROR: Wrong SELECT_CHROM value ($chr). Some values are greater than CHROM_SET variable."
      COUNT_ERROR=$((COUNT_ERROR + 1))
      echo
   fi
done
wait


### INCLUDE_CHROM_X
# Defaults
if [[ -z "$INCLUDE_CHROM_X" ]]; then 
   export INCLUDE_CHROM_X=no
fi

if [[ "${INCLUDE_CHROM_X,,}" != "no" && "${INCLUDE_CHROM_X,,}" != "yes" ]]; then
   echo "#@ WARNING: Invalid option in parameter file (INCLUDE_CHROM_X variable) "
   echo "#@ IVDP will consider INCLUDE_CHROM_X=no"
   export INCLUDE_CHROM_X=no
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
else
   export INCLUDE_CHROM_X=${INCLUDE_CHROM_X,,}
fi


### INCLUDE_CHROM_Y
# Defaults
if [[ -z "$INCLUDE_CHROM_Y" ]]; then 
   export INCLUDE_CHROM_Y=no
fi

if [[ "${INCLUDE_CHROM_Y,,}" != "no" && "${INCLUDE_CHROM_Y,,}" != "yes" ]]; then
   echo "#@ WARNING: Invalid option in parameter file (INCLUDE_CHROM_Y variable) "
   echo "#@ IVDP will consider INCLUDE_CHROM_Y=no"
   export INCLUDE_CHROM_Y=no
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
else
   export INCLUDE_CHROM_Y=${INCLUDE_CHROM_Y,,}
fi


### INCLUDE_CHROM_MT
# Defaults
if [[ -z "$INCLUDE_CHROM_MT" ]]; then 
   export INCLUDE_CHROM_MT=no
fi

if [[ "${INCLUDE_CHROM_MT,,}" != "no" && "${INCLUDE_CHROM_MT,,}" != "yes" ]]; then
   echo "#@ WARNING: Invalid option in parameter file (INCLUDE_CHROM_MT variable) "
   echo "#@ IVDP will consider INCLUDE_CHROM_MT=no"
   export INCLUDE_CHROM_MT=no
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
else
   export INCLUDE_CHROM_MT=${INCLUDE_CHROM_MT,,}
fi


### ADAPTER TRIMMING
if [[ -z "$ADAPTER_TRIMMOMATIC" ]]; then
   #echo "WARNING: IVDP will use the standard TRIMMOMATIC adapter file in $IVDP_DIR/program/03.qc_trim folder"
   #echo "To use a personalized adapter file, specify in parameter file (ADAPTER_TRIMMOMATIC variable)
   #COUNT_WARNING=$((COUNT_WARNING + 1))
   #echo
   if [[ "$TYPE" == "se" ]]; then
      export ADAPTER_TRIMMOMATIC="$IVDP_DIR"/program/03.qc_trim/adapters/TruSeq3-SE.fa:3:30:10
   else
      export ADAPTER_TRIMMOMATIC="$IVDP_DIR"/program/03.qc_trim/adapters/TruSeq3-PE.fa:2:30:10:8:true
   fi
else
   ADAPTER_TRIMMOMATIC=$(ADAPTER_TRIMMOMATIC,,)
fi


### MIN_READ_LENGTH
if [[ -z "$MIN_READ_LENGTH" ]]; then
   export MIN_READ_LENGTH=36
else
   export MIN_READ_LENGTH=$(echo $MIN_READ_LENGTH | tr -d [:alpha:])
fi


### MAX_READ_LENGTH
if [[ -z "$MAX_READ_LENGTH" ]]; then
   export MAX_READ_LENGTH=150
else
   export MAX_READ_LENGTH=$(echo $MAX_READ_LENGTH | tr -d [:alpha:])
fi


### DOWNSAMPLING_FASTQ
if [[ -z "$DOWNSAMPLING_FASTQ" ]]; then
   export DOWNSAMPLING_FASTQ=0
else
   # Check if DOWNSAMPLING_FASTQ is integer or decimal
   # https://stackoverflow.com/questions/14475170/how-do-i-validate-decimal-numbers
   if [[ "$DOWNSAMPLING_FASTQ" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
      export DOWNSAMPLING_FASTQ=$DOWNSAMPLING_FASTQ
   else
      echo "#@ ERROR: Use only numeric values for DOWNSAMPLING_FASTQ variable in parameter file "
      echo "#@ DOWNSAMPLING_FASTQ=0 means no reads sampling from fastq files."
      COUNT_ERROR=$((COUNT_ERROR + 1))
      echo
   fi
fi


### DOWNSAMPLING_BAM
if [[ -z "$DOWNSAMPLING_BAM" ]]; then
   export DOWNSAMPLING_BAM=0
else
   # Check if DOWNSAMPLING_BAM is integer or decimal
   # https://stackoverflow.com/questions/14475170/how-do-i-validate-decimal-numbers
   if [[ "$DOWNSAMPLING_BAM" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
      export DOWNSAMPLING_BAM=$DOWNSAMPLING_BAM
   else
      echo "#@ ERROR: Use only numeric values for DOWNSAMPLING_BAM variable in parameter file "
      echo "#@ DOWNSAMPLING_BAM=0 means no reads sampling from fastq files."
      COUNT_ERROR=$((COUNT_ERROR + 1))
      echo
   fi
fi


### DOWNSAMPLING_FASTQ AND DOWNSAMPLING_BAM
if [[ "$(echo "$DOWNSAMPLING_FASTQ > 0" | bc -l)" -gt 0 && "$(echo "$DOWNSAMPLING_BAM > 0" | bc -l)" -gt 0 ]]; then
   echo "#@ ERROR: use DOWNSAMPLING_FASTQ or DOWNSAMPLING_BAM. Both options cannot be activate at the same time"
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi


### MIN_DEPTH
if [[ "$STEP_VCF_FILTER" == "yes" ]]; then
   if [[ -z "$MIN_DEPTH" ]]; then
      export MIN_DEPTH=3
   else
      export MIN_DEPTH=$(echo $MIN_DEPTH | tr -d [:alpha:])
   fi
else
   export MIN_DEPTH=0
fi


### MAF
if [[ "$STEP_VCF_FILTER" == "yes" ]];then
   if [[ -z "$MAF" ]]; then
      export MAF=0.005
   else
      export MAF=$(echo $MAF | tr -d [:alpha:])
   fi
else
   export MAF=0
fi      


### MISSING
if [[ "$STEP_VCF_FILTER" == "yes" ]];then
   if [[ -z "$MISSING" ]]; then
      export MISSING=0.3
   else
      export MISSING=$(echo $MISSING | tr -d [:alpha:])
   fi
else
   export MISSING=1
fi     


### COMBINE VCF FILES
# Defaults
if [[ "$STEP_VCF_FILTER" == "yes" ]]; then

   COUNT_PROGRAMS=0
   for i in $CALLER_PROGRAM; do COUNT_PROGRAMS=$((COUNT_PROGRAMS + 1)); done
   
   if [[ "$COUNT_PROGRAMS" -gt 1 ]]; then
      
      if [[ -z "$COMBINE_VCF" ]]; then
         export COMBINE_VCF=partial,full
      else
         COMBINE_VCF_SET="none partial full"
         COMBINE_VCF=${COMBINE_VCF,,}
         
         COUNT_COMBINE=0
         COUNT_MATCH=0
         for i in $COMBINE_VCF; do
            COUNT_COMBINE=$((COUNT_COMBINE + 1))
         
            for j in $COMBINE_VCF_SET; do
               if [[ $i == $j ]]; then
                  COUNT_MATCH=$((COUNT_MATCH + 1))
               fi
         
            done
         done
         
         if [[ "$COUNT_COMBINE" != "$COUNT_MATCH" ]]; then
            echo "#@ ERROR: Invalid options for COMBINE_VCF variable in parameter file "
            echo "#@ if STEP_VCF_FILTER=yes, choose COMBINE_VCF=none or COMBINE_VCF=partial or COMBINE_VCF=full or COMBINE_VCF=partial, full"
            COUNT_ERROR=$((COUNT_ERROR + 1))
         else
            export COMBINE_VCF=${COMBINE_VCF,,}
         fi
      fi
   fi
   
else
   export COMBINE_VCF=none
fi


### NUMBER OF THREADS
#https://stackoverflow.com/questions/394230/how-to-detect-the-os-from-a-bash-script
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
   if [[ -z "$THREADS" ]]; then
      NPROC=$(nproc)
      if [[ "$NPROC" -le 16 ]]; then
         export THREADS=$((NPROC - 1))
      else
         export THREADS=16
      fi
      
   else
      export THREADS=$(echo $THREADS | tr -d [:alpha:])
   fi

elif [[ "$OSTYPE" == "darwin"* ]]; then
   if [[ -z "$THREADS" ]]; then
      NPROC=$(sysctl -n hw.ncpu)
      if [[ "$NPROC" -le 16 ]]; then
         export THREADS=$((NPROC - 1))
      else
         export THREADS=16
      fi
      
   else
      export THREADS=$(echo $THREADS | tr -d [:alpha:])
   fi

else
   echo "#@ ERROR: This is not a linux or mac machine. The analysis will stop "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo

fi


### NUMBER OF BATCH SAMPLES
if [[ -z "$BATCH" ]]; then
   export BATCH=5
else
   if [[ "$ANALYSIS" == "1" && "$(echo $BATCH | tr -d [:alpha:])" -gt "20" ]]; then
      export BATCH=20
      
   elif [[ "$ANALYSIS" == "2" && "$(echo $BATCH | tr -d [:alpha:])" -gt "20" ]]; then
      export BATCH=20
   else
      export BATCH=$(echo $BATCH | tr -d [:alpha:])
   fi
fi


### KEEP LOG FILES
# Defaults
if [[ -z "$KEEP_LOG" ]]; then 
   export KEEP_LOG=yes
fi

if [[ "${KEEP_LOG,,}" != "no" && "${KEEP_LOG,,}" != "yes" ]]; then
   echo "#@ WARNING: Invalid option in parameter file variable (KEEP_LOG variable) "
   echo "#@ IVDP will consider KEEP_LOG=no"
   export KEEP_LOG=yes
   COUNT_WARNING=$((COUNTparameter_WARNING + 1))
   echo
else
   export KEEP_LOG=${KEEP_LOG,,}
fi


### KEEP INTERMEDIATE DATA
# Defaults
if [[ -z "$KEEP_INTERMEDIATE_DATA" ]]; then 
   export KEEP_INTERMEDIATE_DATA=yes
fi

if [[ "${KEEP_INTERMEDIATE_DATA,,}" != "no" && "${KEEP_INTERMEDIATE_DATA,,}" != "yes" ]]; then
   echo "#@ WARNING: Invalid option in parameter file variable (KEEP_INTERMEDIATE_DATA variable) "
   echo "#@ IVDP will consider KEEP_INTERMEDIATE_DATA=no"
   export KEEP_INTERMEDIATE_DATA=yes
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
else
   export KEEP_INTERMEDIATE_DATA=${KEEP_INTERMEDIATE_DATA,,}
fi

###############################################
########## CHECK ERRORS AND WARNINGS ##########
###############################################

echo "***** #@ TOTAL OF ERRORS: $COUNT_ERROR *****"
echo "***** #@ TOTAL OF WARNINGS: $COUNT_WARNING *****"

if [[ "$COUNT_ERROR" -gt "0" ]]; then
   echo
   echo "#@ Solve all the ERRORS before run IVDP again "
   echo "Don't worry about WARNINGS: IVDP will use default parameters "
   echo "Aborting analysis"
   exit 1
else
   if [[ "$COUNT_WARNING" -gt "0" ]]; then
      echo
      echo "Don't worry about WARNINGS: IVDP will use default parameters "
   fi
fi

#echo
#echo "#@#############################################################"

