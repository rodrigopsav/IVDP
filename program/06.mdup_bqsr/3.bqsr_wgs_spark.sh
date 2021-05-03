#!/bin/bash


echo "#@##################################################################"
echo "#@ BASE QUALITY SCORE RECALIBRATION SPARK: SAMPLE ${i}"
echo
D1=`date "+%D    %T"`
echo "#@ Date and Time: $D1"
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

gatk --java-options '-Xmx20g -XX:+UseParallelGC -XX:ParallelGCThreads=4' BaseRecalibratorSpark \
	-R $REFERENCE \
	-I $ALIGN_MDUP/${i}_trim_"$ALIGNER"_sorted_RG_mdup.bam \
	-known-sites $VARIANTS \
	-O $ALIGN_BQSR/${i}_recalibration.table \
	--tmp-dir $TMP \
	-- --spark-runner LOCAL --spark-master local[$THREADS_SPARK] \
	--conf spark.local.dir=$TMP
	

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ BQSR-SPARK SAMPLE ${i} TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "#@##################################################################"

echo
echo
echo

echo "#@##################################################################"
echo "#@ APPLY BQSR SPARK: SAMPLE ${i}"
echo
D2=`date "+%D    %T"`
echo "#@ Date and Time: $D2"
echo
       #################################


START2=$(date +%s)

gatk --java-options '-Xmx20g -XX:+UseParallelGC -XX:ParallelGCThreads=4' ApplyBQSRSpark \
	-R $REFERENCE \
	-I $ALIGN_MDUP/${i}_trim_"$ALIGNER"_sorted_RG_mdup.bam \
	--bqsr-recal-file $ALIGN_BQSR/${i}_recalibration.table \
	--static-quantized-quals 10 --static-quantized-quals 20 \
	--static-quantized-quals 30 \
	-O $ALIGN_BQSR/${i}_trim_"$ALIGNER"_sorted_RG_mdup_bqsr.bam \
	--tmp-dir $TMP \
	-- --spark-runner LOCAL --spark-master local[$THREADS_SPARK] \
	--conf spark.local.dir=$TMP
wait

rm -r $ALIGN_BQSR/${i}_recalibration.table 2> /dev/null
rm -r $ALIGN_BQSR/${i}_trim_"$ALIGNER"_sorted_RG_mdup_bqsr.bam.sbi 2> /dev/null
wait

END2=$(date +%s)
DIFF2=$(( $END2 - $START2 ))

echo
echo "#@ APPLY BQSR-SPARK SAMPLE ${i} TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF2/3600)) $(($DIFF2%3600/60)) $(($DIFF2%60)))"
echo
echo "#@##################################################################"

echo
