#!/bin/bash


echo "#@#############################################################"
echo "#@ ALIGNMENT STATISTICS: $OUTPUT_NAME"
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo

START1=$(date +%s)

### CHECK IF INITIAL FILES EXIST

if [ ! -f "$REPORT_ALIGN/"$ALIGNER"_align_stats.txt" ]; then
   touch $REPORT_ALIGN/"$ALIGNER"_align_stats.txt
   echo $'Sample\tDepth_Cov\tBreadth_Cov\tDepth_Only_Cov_Regions\tTotal_Reads\tReads_Mapped\tReads_Unmapped\tReads_Paired\tBases_Mapped\tRead_Length\tBase_Quality\tMapped_Reads_Perc\tUnmapped_Reads_Perc' > $REPORT_ALIGN/"$ALIGNER"_align_stats.txt
fi
wait


#if [ ! -f "$REPORT_ALIGN/"$ALIGNER"_depth_by_chrom.txt" ]; then
#   touch $REPORT_ALIGN/"$ALIGNER"_depth_by_chrom.txt
#   paste <(echo "Samples") <(echo "$(awk 'BEGIN { ORS="\t" } { print }' $CHROM_NAME)") --delimiters '\t' > $REPORT_ALIGN/"$ALIGNER"_depth_by_chrom.txt
#fi
#wait


#if [ ! -f "$REPORT_ALIGN/"$ALIGNER"_breadth_by_chrom.txt" ]; then
#   touch $REPORT_ALIGN/"$ALIGNER"_breadth_by_chrom.txt
#   paste <(echo "Samples") <(echo "$(awk 'BEGIN { ORS="\t" } { print }' $CHROM_NAME)") --delimiters '\t' > $REPORT_ALIGN/"$ALIGNER"_breadth_by_chrom.txt
#fi
#wait


###############################
##### PREPARE STATS FILES #####
###############################

### SELECT SPECIFIC LINES OF STAT FILE
#for lines in 1 7 9 11 20 25 31; do
#   cat $STATS_ALIGN/"${i}".stats | grep ^SN | cut -f 2- | \
#   awk -v var="$lines" -F"\t" 'NR==var {print $2}' >> $STATS_ALIGN/"${i}".tmp1
#done
#wait
grep ^SN $STATS_ALIGN/"${i}".stats | grep "raw total sequences:" | cut -f 3 >> $STATS_ALIGN/"${i}".tmp1
grep ^SN $STATS_ALIGN/"${i}".stats | grep "reads mapped:" | cut -f 3 >> $STATS_ALIGN/"${i}".tmp1
grep ^SN $STATS_ALIGN/"${i}".stats | grep "reads unmapped:" | cut -f 3 >> $STATS_ALIGN/"${i}".tmp1
grep ^SN $STATS_ALIGN/"${i}".stats | grep "reads paired:" | cut -f 3 >> $STATS_ALIGN/"${i}".tmp1
grep ^SN $STATS_ALIGN/"${i}".stats | grep "bases mapped (cigar):" | cut -f 3 >> $STATS_ALIGN/"${i}".tmp1
grep ^SN $STATS_ALIGN/"${i}".stats | grep "average length:" | cut -f 3 >> $STATS_ALIGN/"${i}".tmp1
grep ^SN $STATS_ALIGN/"${i}".stats | grep "average quality:" | cut -f 3 >> $STATS_ALIGN/"${i}".tmp1



### TRANSPOSE 
awk 'BEGIN { ORS="\t" } { print }' $STATS_ALIGN/"${i}".tmp1 > $STATS_ALIGN/"${i}".tmp2
wait


### CREATE 2 NEW COLUMNS: PERCENTAGE OF MAPPED AND UNMAPPED READS
awk -v OFS='\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\n", $1, $2, $3, $4, $5, $6, $7, $2/$1*100, $3/$1*100}' $STATS_ALIGN/${i}.tmp2 > $STATS_ALIGN/${i}.tmp2."$$" && mv $STATS_ALIGN/${i}.tmp2."$$" $STATS_ALIGN/${i}.tmp2
wait


### PASTE SAMPLES AND TRANSPOSED FILE
if ! grep -q $i $REPORT_ALIGN/"$ALIGNER"_align_stats.txt; then
   paste <(echo "$i") <(echo "$(cat $STATS_ALIGN/"${i}".cov)") <(echo "$(cat $STATS_ALIGN/"${i}".tmp2)") --delimiters '\t' >> $REPORT_ALIGN/"$ALIGNER"_align_stats.txt
fi
wait


### DELETE TEMP FILES
rm $STATS_ALIGN/"${i}".stats
rm $STATS_ALIGN/"${i}".cov
rm $STATS_ALIGN/"${i}".tmp1
rm $STATS_ALIGN/"${i}".tmp2
wait


###############################
##### PREPARE STATS FILES #####
#####       BY CHROM      #####
###############################


### DEPTH AND BREADTH OF COVERAGE
#for chr in $(cat $CHROM_NAME); do
#   SEARCH=$(echo "SN:$chr")
#   totalBasesGenome=$(samtools view -H $ALIGN/"${i}"_trim_"$ALIGNER"_sorted_RG.bam | grep -P '^@SQ' | grep -w $SEARCH | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}')
#   awk -v chr=$chr '$1 == chr {numberPositions = $3 - $2; basesCovered = numberPositions * $4; print numberPositions, basesCovered}' $STATS_ALIGN/${i}.coverageBed | awk -v totalBasesGenome=$totalBasesGenome '$2 > 0 {sumCoverage+=$2; totalBasesCovered+=$1} END {printf "%.2f\t%.2f\t%.2f\n", sumCoverage/totalBasesGenome, (totalBasesCovered/totalBasesGenome)*100, sumCoverage/totalBasesCovered}' >> $STATS_ALIGN/"${i}".tmp3
#done
#wait

#awk -F"\t" '{print $1}' $STATS_ALIGN/"${i}".tmp3 | awk 'BEGIN { ORS="\t" } { print }' > $STATS_ALIGN/"${i}".tmp4
#wait
#awk -F"\t" '{print $2}' $STATS_ALIGN/"${i}".tmp3 | awk 'BEGIN { ORS="\t" } { print }' > $STATS_ALIGN/"${i}".tmp5
#wait

#paste <(echo "$i") <(echo "$(cat $STATS_ALIGN/"${i}".tmp4)") --delimiters '\t' >> $REPORT_ALIGN/"$ALIGNER"_depth_by_chrom.txt
#wait
#paste <(echo "$i") <(echo "$(cat $STATS_ALIGN/"${i}".tmp5)") --delimiters '\t' >> $REPORT_ALIGN/"$ALIGNER"_breadth_by_chrom.txt
#wait

#rm $STATS_ALIGN/"${i}".tmp3
#rm $STATS_ALIGN/"${i}".tmp4
#rm $STATS_ALIGN/"${i}".tmp5
#wait

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ STATISTICS OF ALIGNMENT TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@ Complete process"
echo "@#############################################################"
