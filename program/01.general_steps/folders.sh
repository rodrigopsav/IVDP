#!/bin/bash

##### TEMPORARY FILES FOLDER
export TMP="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/tmp_"${ANALYSIS_ID}"


##### LOG FOLDER INDEX REFERENCE GENOME
export LOG_INDEX="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/log/index_reference_genome


##### FOLDERS FOR DATA
#export SRA_DIR="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/sra_fastq
export TRIM="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/fastqQC
export ALIGN="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/alignment/raw/"$ALIGNER"
export ALIGN_MDUP="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/alignment/mdup/"$ALIGNER"
export ALIGN_BQSR="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/alignment/bqsr/"$ALIGNER"
export GENECOUNT_STAR="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/gene_count/STAR
export GENECOUNT_DATA="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/gene_count/"$GENE_COUNT"_"$ALIGNER"
export GVCF_RAW="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/gvcf
export VCF_RAW="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/vcf/raw
export VCF_FILTERED_SNP="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/vcf/filtered/snp
export VCF_FILTERED_INDEL="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/vcf/filtered/indel
export VCF_FILTERED_COMBINED="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/vcf/combined


#### FOLDERS FOR LOGS
export LOG_TRIM="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/log/fastqQC
export LOG_ALIGN="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/log/alignment/raw/"$ALIGNER"
export LOG_ALIGN_MDUP="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/log/alignment/mdup/"$ALIGNER"
export LOG_ALIGN_BQSR="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/log/alignment/bqsr/"$ALIGNER"
export LOG_GENECOUNT="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/log/gene_count
export LOG_VCF_RAW="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/log/vcf/raw
export LOG_VCF_FILTERED="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/log/vcf/filtered


##### FOLDERS FOR STATS
export STATS_QC1="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/stats/fastqQC/fastq_raw
export STATS_QC2="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/stats/fastqQC/fastq_trimmed
export STATS_TRIM="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/stats/fastqQC/trimmomatic
export STATS_ALIGN="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/stats/alignment/"$ALIGNER"
export STATS_VCF_RAW="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/stats/vcf/raw
export STATS_VCF_FILTERED_SNP="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/stats/vcf/filtered/snp
export STATS_VCF_FILTERED_INDEL="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/stats/vcf/filtered/indel
export STATS_VCF_FILTERED_COMBINED="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/stats/vcf/combined


##### FOLDERS FOR MULTIQC REPORT
export REPORT_TRIM="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/report/fastqQC
export REPORT_ALIGN="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/report/alignment
export REPORT_GENECOUNT="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/report/gene_count
export REPORT_VCF_RAW="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/report/vcf/raw
export REPORT_VCF_FILTERED_SNP="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/report/vcf/filtered/snp
export REPORT_VCF_FILTERED_INDEL="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/report/vcf/filtered/indel
export REPORT_VCF_FILTERED_COMBINED="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/report/vcf/combined


######################################
