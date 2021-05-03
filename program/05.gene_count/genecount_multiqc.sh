#!/bin/bash


echo "#@##################################################################"
echo "#@ GENE COUNT HTSEQ REPORT"
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo
       #################################


START1=$(date +%s)

multiqc -f -n "$FEATURE"_"$GENE_COUNT"_"$ALIGNER".html $GENECOUNT_DATA/*_counts_"$FEATURE".txt -o $REPORT_GENECOUNT

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ GENE COUNT HTSEQ STATISTICS TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "#@##################################################################"

echo
echo
echo

