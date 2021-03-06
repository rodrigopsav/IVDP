#!/bin/bash

TRIM="$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/data/fastqQC
if [[ $RAW_DATA == $TRIM ]]; then
   R1=${i}.R1_paired.fastq.gz
   R2=${i}.R2_paired.fastq.gz

else
   R1=${i}${EXTENSION_PE1}
   R2=${i}${EXTENSION_PE2}
fi
wait


echo "#@##################################################################"
echo "#@ SAMPLE: $i --- MAPPING $ALIGNER"
echo
D1=`date "+%D    %T"`
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

cd $RAW_DATA
bbmap.sh in=$R1 in2=$R2 out=$ALIGN/${i}_trim_"$ALIGNER".sam path=$bbmap_index threads=$THREADS maxindel=16000
wait 

#> $ALIGN/${i}_trim_"$ALIGNER".sam
#wait
	
### SAM TO BAM
sambamba view -S -f bam -t $THREADS $ALIGN/${i}_trim_"$ALIGNER".sam > $ALIGN/${i}_trim_"$ALIGNER".bam
wait

rm -r $ALIGN/${i}_trim_"$ALIGNER".sam

### SORTED BAM

#https://stackoverflow.com/questions/394230/how-to-detect-the-os-from-a-bash-script
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
   if [[ $(nproc) -le 8 ]]; then
      export NPROC=$(nproc)
   else
      export NPROC=8
   fi
elif [[ "$OSTYPE" == "darwin"* ]]; then
   if [[ $(sysctl -n hw.ncpu) -le 8 ]]; then
      export NPROC=$(sysctl -n hw.ncpu)
   else
      export NPROC=8
   fi
else
   echo""
fi

sambamba sort -t $NPROC --tmpdir="$TMP" -o $ALIGN/${i}_trim_"$ALIGNER"_sorted.bam $ALIGN/${i}_trim_"$ALIGNER".bam
wait

rm -r $ALIGN/${i}_trim_"$ALIGNER".bam

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ MAPPING SAMPLE $i TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "#@#############################################################"

echo
echo
echo
