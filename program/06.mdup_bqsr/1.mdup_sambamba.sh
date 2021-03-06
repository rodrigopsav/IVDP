#!/bin/bash


echo "#@##################################################################"
echo "#@ MARK DUPLICATED READS: SAMPLE $i"
echo
D1=`date "+%D    %T"`
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

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

sambamba markdup -t $NPROC --overflow-list-size 1000000 \
--hash-table-size 1000000 --tmpdir="$TMP" \
$ALIGN/${i}_trim_"$ALIGNER"_sorted_RG.bam \
$ALIGN_MDUP/${i}_trim_"$ALIGNER"_sorted_RG_mdup.bam
wait

#rm -v $ALIGN/${i}_trim_"$ALIGNER"_sorted_RG.bam
#rm -v $ALIGN/${i}_trim_"$ALIGNER"_sorted_RG.bam.bai


END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ MDUP SAMPLE $i TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "#@##################################################################"

echo

