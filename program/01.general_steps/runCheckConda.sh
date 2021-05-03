#!/bin/bash

echo
echo
echo "#@############################################################"
echo "#@                      CHECK CONDA ENVIRONMENTS"
#echo


#############################################
########## LOAD DEFAULT PARAMETERS ##########
#############################################
COUNT_ERROR=0
COUNT_WARNING=0

###
if ! command -v parallel &> /dev/null; then
   echo "#@ ERROR: parallel cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v dos2unix &> /dev/null; then
   echo "#@ ERROR: dos2unix cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v rename &> /dev/null; then
   echo "#@ ERROR: rename cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v sed &> /dev/null; then
   echo "#@ ERROR: sed cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v git &> /dev/null; then
   echo "#@ ERROR: git cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v pigz &> /dev/null; then
   echo "#@ ERROR: pigz cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v datamash &> /dev/null; then
   echo "#@ ERROR: datamash cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v prefetch &> /dev/null; then
   echo "#@ ERROR: sra-toolkit cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v esearch &> /dev/null; then
   echo "#@ ERROR: entrez-direct cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v fastqc &> /dev/null; then
   echo "#@ ERROR: fastqc cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v fastq_info &> /dev/null; then
   echo "#@ ERROR: fastq_utils cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v seqtk &> /dev/null; then
   echo "#@ ERROR: seqtk cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v multiqc &> /dev/null; then
   echo "#@ ERROR: multiqc cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v trimmomatic &> /dev/null; then
   echo "#@ ERROR: trimmomatic cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v bbmap.sh &> /dev/null; then
   echo "#@ ERROR: bbmap.sh cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v bowtie &> /dev/null; then
   echo "#@ ERROR: bowtie cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v bowtie2 &> /dev/null; then
   echo "#@ ERROR: bowtie2 cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v bwa &> /dev/null; then
   echo "#@ ERROR: bwa cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v bwa-mem2 &> /dev/null; then
   echo "#@ ERROR: bwa-mem2 cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v gsnap &> /dev/null; then
   echo "#@ ERROR: gsnap cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v hisat2 &> /dev/null; then
   echo "#@ ERROR: hisat2 cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v STAR &> /dev/null; then
   echo "#@ ERROR: STAR cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v sambamba &> /dev/null; then
   echo "#@ ERROR: sambamba cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v samtools &> /dev/null; then
   echo "#@ ERROR: samtools cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v bgzip &> /dev/null; then
   echo "#@ ERROR: htslib cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v picard &> /dev/null; then
   echo "#@ ERROR: picard cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v bcftools &> /dev/null; then
   echo "#@ ERROR: bcftools cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v gatk &> /dev/null; then
   echo "#@ ERROR: gatk cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v varscan &> /dev/null; then
   echo "#@ ERROR: varscan cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v freebayes &> /dev/null; then
   echo "#@ ERROR: freebayes cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v freebayes-parallel &> /dev/null; then
   echo "#@ ERROR: freebayes-parallel cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v platypus &> /dev/null; then
   echo "#@ ERROR: platypus cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v vcftools &> /dev/null; then
   echo "#@ ERROR: vcftools cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v plink &> /dev/null; then
   echo "#@ ERROR: plink cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v plink2 &> /dev/null; then
   echo "#@ ERROR: plink2 cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v bedtools &> /dev/null; then
   echo "#@ ERROR: bedtools cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v htseq-count &> /dev/null; then
   echo "#@ ERROR: htseq-count cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v featureCounts &> /dev/null; then
   echo "#@ ERROR: featureCounts cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v snpEff &> /dev/null; then
   echo "#@ ERROR: snpEff cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v SnpSift &> /dev/null; then
   echo "#@ ERROR: SnpSift cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v vcffilter &> /dev/null; then
   echo "#@ ERROR: vcflib cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v gs &> /dev/null; then
   echo "#@ ERROR: ghostscript cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi

###
if ! command -v R &> /dev/null; then
   echo "#@ ERROR: R Program cannot be found. "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi


##################################
########## CHECK ERRORS ##########
##################################
echo "***** #@ TOTAL OF ERRORS: $COUNT_ERROR *****"

if [[ "$COUNT_ERROR" -gt "0" ]]; then
   echo
   echo "#@ Install IVDP dependencies before run IVDP again "
   exit 1
fi

#echo
#echo "#@#############################################################"


