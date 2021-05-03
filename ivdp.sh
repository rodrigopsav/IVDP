#!/bin/bash

usage() { echo "Usage: $0 -p params.txt -c configSlurm.txt
       -p: ivdp parameter file
       -c: slurm parameter file (required only for the HPCC with slurm)

Variables in the parameter file

INPUT_DIR       # Directory path for the original fastq files.
OUTPUT_DIR      # Directory path for ivdp output (Default: home directory if OUTPUT_DIR is ommited).
OUTPUT_NAME     # Name for ivdp analysis. Choose a single name without spaces (Default: user name if OUTPUT_NAME is ommited).
LIST_SAMPLES    # File path for a list with samples id and individual id, separated by space. Use LIST_SAMPLES=none if no list is available.
ANALYSIS        # ANALYSIS=1 (for WGS) or ANALYSIS=2 (for RNAseq).
EXTENSION_SE    # For single-end fastq files, write the extension of original fastq files. If paired-end files, left it blank (EXTENSION_SE= ) or remove EXTENSION_SE from the parameters file.
EXTENSION_PE1   # For paired-end fastq files, write the extension of PE1 fastq files. If single-end files, left it blank (EXTENSION_PE1= ) or remove EXTENSION_PE1 from the parameters file.
EXTENSION_PE2   # For paired-end fastq files, write the extension of PE2 fastq files. If single-end files, left it blank (EXTENSION_PE2= ) or remove EXTENSION_PE2 from the parameters file.


ASSEMBLY      # Give a name for the reference genome. Eg ASSEMBLY=hg38.
REFERENCE     # File path for reference genome (*.fa). The reference genome must be uncompressed file and not (*.fa.gz).
VARIANTS      # File path for variants file (*.vcf.gz) for the reference genome. Use VARIANTS=none if no variants file is available. The variants must be (*.vcf.gz) and not (*.vcf).
ANNOTATION    # File path for annotation file (*.gtf) for the reference genome. Use ANNOTATION=none if no annotation file is available. The annotation must be uncompressed file and not (*.gtf.gz).


STEP_QC_TRIM        # Activate quality control and adapter trimming step? STEP_QC_TRIM=yes or STEP_QC_TRIM=no (Default yes)
STEP_ALIGNMENT      # Activate alignment step? STEP_ALIGNMENT=yes or STEP_ALIGNMENT=no (Default yes)
STEP_GENECOUNT      # Activate gene count step? STEP_GENECOUNT=yes or STEP_GENECOUNT=no (Default: yes for RNAseq and no for WGS)
STEP_MDUP           # Activate mark duplicated read step? STEP_MDUP=yes or STEP_MDUP=no (Default: yes). It includes splitncigar step for RNAseq
STEP_BQSR           # Activate base quality score recalibration step? STEP_BQSR=yes or STEP_BQSR=no (Default: yes)
STEP_VAR_CALL       # Activate variant calling step? STEP_VAR_CALL=yes or STEP_VAR_CALL=no (Default: yes)
STEP_GVCF_TO_VCF    # Activate merge gvcfs and convert gvcfs to vcf file (only for CALLER_PROGRAM=gatk4GenomicsDBImport or CALLER_PROGRAM=gatk4CombineGVCFs). STEP_GVCF_TO_VCF=yes or STEP_GVCF_TO_VCF=no (Default: yes)
STEP_VCF_FILTER     # Activate vcf filtering step? STEP_VCF_FILTER=yes or STEP_VCF_FILTER=no (Default: yes)


ALIGNER_PROGRAM     # Choose programs for alignment step. If more than one program are chosen, use ',' between them (Options: bbmap, bowtie, bowtie2, bwa, bwa2, gsnap, hisat2, star). Default: ALIGNER_PROGRAM=bwa2 for WGS (if ANALYSIS=1) and ALIGNER_PROGRAM=hisat2 for RNAseq (if ANALYSIS=2).
CALLER_PROGRAM      # Choose programs for variant calling step. If more than one program are chosen, use ',' between them (Options: bcftools, freebayes, freebayes-parallel, gatk4, gatk4GenomicsDBImport, gatk4CombineGVCFs, platypus, varscan). Default: CALLER_PROGRAM=gatk4GenomicsDBImport for WGS (if ANALYSIS=1) and CALLER_PROGRAM=freebayes for RNAseq (if ANALYSIS=2)
MDUP_PROGRAM        # Choose the program for mark duplicated reads. (Options: sambamba, markdupspark, or picard). Default: sambamba
BQSR_PROGRAM        # Choose the program for base quality score recalibration. (Options: bqsrspark or bqsr). Default: bqsrspark
GENECOUNT_PROGRAM   # Choose programs for gene count step. If more than one program are chosen, use ',' between them (Options: htseq, featurecounts). Use GENECOUNT_PROGRAM=none if you do not want gene counts. Default: GENECOUNT_PROGRAM=none for WGS (if ANALYSIS=1) and GENECOUNT_PROGRAM=htseq for RNAseq (if ANALYSIS=2)
FEATURE_TYPE        # Any feature of the third column of gtf file. (Default FEATURE_TYPE=exon). Use FEATURE_TYPE=none if you do not want gene counts.


BP_BY_CHROM       # Split chromosomes in chunks of bp for variant calling (It is activated only if CALL_BY_CHROM=yes). (Default: BP_BY_CHROM=all (use all bp of the chromosome))
CHROM_SET         # Number of autosome pairs (Eg: CHROM_SET=22 for humans, CHROM_SET=29 for cow, CHROM_SET=18 for pig)
SELECT_CHROM      # Comma-separated list of chromosomes that will be used for variant calling (Eg. SELECT_CHROM=1,2,3). If omitted, IVDP will consider all the autosomes. 
INCLUDE_CHROM_X   # Do the variants call for chromosome X? INCLUDE_CHROM_X=no or INCLUDE_CHROM_X=yes (Default: INCLUDE_CHROM_X=no)
INCLUDE_CHROM_Y   # Do the variants call for chromosome Y? INCLUDE_CHROM_Y=no or INCLUDE_CHROM_Y=yes (Default: INCLUDE_CHROM_Y=no)
INCLUDE_CHROM_MT  # Do the variants call for mitochondrial chromosome? INCLUDE_CHROM_MT=no or INCLUDE_CHROM_MT=yes (Default: INCLUDE_CHROM_MT=no)
MIN_READ_LENGTH   # Minimum read length to be kept after adapter trimming (Default 36)
MAX_READ_LENGTH   # Maximum read length to be kept after adapter trimming (Default 150)
DOWNSAMPLING_FASTQ   # Down sampling the reads of fastq files according the desire coverage. (Default 0). DOWNSAMPLING_FASTQ=0 means no read sampling from fastq files.
DOWNSAMPLING_BAM     # Down sampling the reads of bam files according the desire coverage. (Default 0). DOWNSAMPLING_BAM=0 means no read sampling from bam files.
MIN_DEPTH       # Minimum locus coverage for the filtering step (Default: MIN_DEPTH=3 if STEP_VCF_FILTER=yes and MIN_DEPTH=0 if STEP_VCF_FILTER=no)
MAF             # Minimum allele frequency for the filtering step (Default: MAF=0.001 if STEP_VCF_FILTER=yes and MAF=0 if STEP_VCF_FILTER=no)
MISSING         # Maximum  missing genotype allowed. MISSING=0 means no missing genotypes and MISSING=1 means to allow all missing genotypes (no filter). MISSING=0.1 will exclude loci with more than 10% missing genotypes. (Default: MISSING=0.3 if STEP_VCF_FILTER=yes and MISSING=1 if STEP_VCF_FILTER=no)
COMBINE_VCF     # It works only if STEP_VCF_FILTER=yes and if there is more than one option for CALLER_PROGRAM. Choose:
                # COMBINE_VCF=none : deactivate combine different vcfs;
                # COMBINE_VCF=partial : If it was used more than one program for variant calling, combine common SNPs that appeared AT LEAST in two different vcf files.
                # COMBINE_VCF=full : If it was used more than one program for variant calling, combine common SNPs  that appeared IN ALL the vcf files.
                # COMBINE_VCF=partial, full: It outputs both methods.
STATS_ALIGNMENT # Calculate statistics from alignment step? STATS_ALIGNMENT=yes or STATS_ALIGNMENT=no (Default yes)
KEEP_LOG        # Keep the log files in the output folder? KEEP_LOG=no or KEEP_LOG=yes (Default yes)
KEEP_INTERMEDIATE_DATA  # Keep intermediate fastq and bam files in the output folder? KEEP_INTERMEDIATE_DATA=no or KEEP_INTERMEDIATE_DATA=yes (Default yes)


*** USE THESE VARIABLES ONLY WITH LOCAL MACHINES (WON'T WORK ON HPCC WITH SLURM)
CALL_BY_CHROM   # Should variant calling be done by chromosome? CALL_BY_CHROM=no or CALL_BY_CHROM=yes (Default yes)
THREADS         # Number of threads
BATCH           # Number of samples to be processed at same time (maximum of 20 jobs at once)

       " 1>&2; exit 1; }

while getopts :p:c: option; do
   case "${option}" in
   p) PARAMETERS=${OPTARG};;
   c) CONFIGSLURM=${OPTARG};;
   *) usage;;
   esac
done
shift $((OPTIND -1))


##### Run ivdpLocal or ivdpSlurm #####
export IVDP_DIR=$(dirname $(readlink -f $0))


#
if [[ -z "$PARAMETERS" ]]; then
   echo "ERROR: -p flag is empty "
   echo "Aborting analysis"
   usage
   exit 1

else
   export PARAMETERS=$(readlink -f $PARAMETERS)
   
   if [[ ! -f "$PARAMETERS" ]]; then
      echo "ERROR: parameter file does not exist "
      echo "Aborting analysis"
      usage
      exit 1
   fi
fi
wait


#
if [[ "$CONFIGSLURM" ]]; then
   export CONFIGSLURM=$(readlink -f $CONFIGSLURM)

   if [[ ! -f "$CONFIGSLURM" ]]; then
      echo "ERROR: slurm parameter file does not exist"
      echo "Aborting analysis"
      usage
      exit 1
   fi
fi
wait


#
if [[ ! -f "$CONFIGSLURM" ]]; then
   source $IVDP_DIR/program/ivdpLocalRun.sh
else
   for par in $(cat $CONFIGSLURM); do export $par; done
   source $IVDP_DIR/program/ivdpSlurmRun.sh
fi
