#!/bin/bash


######################################
##### CHECK PROGRAMS AND FOLDERS #####
######################################

##### Check Output Folders for This Step
##

source "$IVDP_DIR"/program/01.general_steps/folders.sh
mkdir -p $LOG_INDEX 2> /dev/null
mkdir -p $TMP 2> /dev/null

for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER
   export eval "$ALIGNER"_index="$IVDP_DIR"/dataBases/genome_idx/"$ASSEMBLY"_index_"$ALIGNER"
   mkdir -p $(eval echo \$${ALIGNER}_index)
done

#
#
#

echo
echo
echo "#@#############################################################"
echo "#@ INDEX REFERENCE GENOME "
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo


START1=$(date +%s)


### INDEXING GENOME
for ALIGNER in $ALIGNER_PROGRAM; do
   export ALIGNER=$ALIGNER

   if [[ "$(ls $(eval echo \$${ALIGNER}_index) | wc -l)" != 0 ]] && [[ "$(du -s $(eval echo \$${ALIGNER}_index) | awk '{ print $1 }')" -gt "100000" ]]; then
      echo "${ASSEMBLY}_index_${ALIGNER} ALREADY EXIST"
      echo "REFERENCE INDEX WAS CREATED IN $(eval echo \$${ALIGNER}_index) "

   else

      if [[ "${ANNOTATION,,}" == "none" && "$ALIGNER" == "hisat2" ]]; then
         $IVDP_DIR/program/02.genome_index/hisat2_index_without_annotation.sh > $LOG_INDEX/"$ALIGNER"_index.txt 2>&1
      elif [[ "${ANNOTATION,,}" == "none" && "$ALIGNER" == "star" ]]; then
         $IVDP_DIR/program/02.genome_index/star_index_without_annotation.sh > $LOG_INDEX/"$ALIGNER"_index.txt 2>&1
      else
         $IVDP_DIR/program/02.genome_index/"$ALIGNER"_index.sh > $LOG_INDEX/"$ALIGNER"_index.txt 2>&1 
      fi
      wait
      echo "REFERENCE INDEX WAS CREATED IN $(eval echo \$${ALIGNER}_index) "
   fi
done
wait


### INDEXING REFERENCE FILE
if test -f "$(echo $REFERENCE.fai)"; then 
   echo "INDEXED REFERENCE.FAI ALREADY EXISTS"
else
   samtools faidx $REFERENCE
   echo "INDEXED REFERENCE.FAI WAS CREATED SUCCESSFULLY"
fi
wait


### INDEXING VARIANTS FILE
if test -f "$(echo $VARIANTS.tbi)"; then 
   echo "INDEXED VARIANT.TBI ALREADY EXISTS"
else
   bcftools index -f -t --threads $THREADS $VARIANTS
   echo "INDEXED VARIANT.TBI WAS CREATED SUCCESSFULLY"
fi
wait


### CREATING SEQUENCE DICTIONARY
if test -f "$(echo "${REFERENCE%.fa}".dict)"; then 
   echo "SEQUENCE DICTIONARY $ASSEMBLY ALREADY EXIST"
else
   $IVDP_DIR/program/02.genome_index/dictionary.sh
   echo "SEQUENCE DICTIONARY $ASSEMBLY WAS CREATED SUCCESSFULLY"
fi
wait


### PREPARING FILES FOR VCF ANNOTATION
if [[ "$ANNOTATION" != "none" ]]; then

   # Check if snpEff.conf exist
   source $IVDP_DIR/program/01.general_steps/folders.sh
   SNPEFF_CONFIG=$(echo $IVDP_DIR/dataBases/annotation_db/snpEff.config)
   
   if [[ ! -f "$SNPEFF_CONFIG" ]]; then
      export SNPEFF_FILE=$(readlink -f $(which snpEff))
      export SNPEFF_FOLDER=$(echo "${SNPEFF_FILE%/*}")
      mkdir -p $IVDP_DIR/dataBases/annotation_db
      cp -r $SNPEFF_FOLDER/snpEff.config $IVDP_DIR/dataBases/annotation_db/snpEff.config
   fi
   wait
   
   # Build snpEff database
   source $IVDP_DIR/program/01.general_steps/folders.sh
   SNPEFF_DATA=$(echo $IVDP_DIR/dataBases/annotation_db/data/"$ASSEMBLY")
   #if [[ -d "DIR" ]] tests if directory exist
   #if [[ ! -d "DIR" ]] tests if a directory does not exist
   if [[ ! -d "$SNPEFF_DATA" ]]; then
      echo "CREATING $ASSEMBLY SNPEFF DATABASE"
      echo
      $IVDP_DIR/program/08.vcf_annotation/snpEff_build.sh
   fi
   wait
fi
wait


END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ INDEX REFERENCE GENOME TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@#############################################################"

