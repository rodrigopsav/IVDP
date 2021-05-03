#!/bin/bash


echo "#@##################################################################"
echo "#@ MULTIQC RAW FILE"
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo
       #################################


START1=$(date +%s)

multiqc -f -n report_fastq_raw.html $STATS_QC1 -o $REPORT_TRIM

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ MULTIQC RAW TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "#@##################################################################"

echo
echo
echo

