#!/bin/bash


echo "@##################################################################"
echo "@ MARK DUPLICATED READS SPARK: SAMPLE $i"
echo
D1=`date "+%D    %T"`
echo "@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

#https://stackoverflow.com/questions/394230/how-to-detect-the-os-from-a-bash-script
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
   if [[ $(nproc) -le 8 ]]; then
      export THREADS_SPARK=$(nproc)
   else
      export THREADS_SPARK=8
   fi
elif [[ "$OSTYPE" == "darwin"* ]]; then
   if [[ $(sysctl -n hw.ncpu) -le 8 ]]; then
      export THREADS_SPARK=$(sysctl -n hw.ncpu)
   else
      export THREADS_SPARK=8
   fi
else
   echo""
fi

gatk --java-options '-Djava.io.tmpdir=$TMP -Xmx32g' MarkDuplicatesSpark \
   -I $ALIGN/${i}_trim_"$ALIGNER"_sorted_RG.bam \
   -O $ALIGN_MDUP/${i}_trim_"$ALIGNER"_sorted_RG_mdup.bam \
   --tmp-dir $TMP \
   -- --spark-runner LOCAL --spark-master local[$THREADS_SPARK] \
   --conf spark.local.dir=$TMP
wait

#rm -v $ALIGN/${i}_trim_"$ALIGNER"_sorted_RG.bam
#rm -v $ALIGN/${i}_trim_"$ALIGNER"_sorted_RG.bam.bai


END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "@ MDUP SAMPLE $i TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "@##################################################################"

echo

