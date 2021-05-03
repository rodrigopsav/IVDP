#!/bin/bash


echo "#@##################################################################"
echo "#@ ALIGNMENT REPORT"
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

if [[ $ALIGNER == "star" ]]; then
   multiqc -f -n star.html $STATS_ALIGN/*Log.final.out -o $REPORT_ALIGN
   multiqc -f -n exon_star.html $GENECOUNT_STAR/*ReadsPerGene.out.tab -o $REPORT_GENECOUNT
else
   multiqc -f -n "$ALIGNER"_align_multiqc.html $STATS_ALIGN/*.txt -o $REPORT_ALIGN
fi
wait

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ ALIGNMENT STATISTICS TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "#@##################################################################"

echo
echo
echo

