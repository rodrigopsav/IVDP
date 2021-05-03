#!/bin/bash

echo
echo "#@#############################################################"
echo "#@ STAR1PASS-INDEX $ASSEMBLY"
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

# https://software.broadinstitute.org/gatk/documentation/article.php?id=3891

STAR \
--runMode genomeGenerate \
--runThreadN $THREADS \
--genomeDir $star_index \
--genomeFastaFiles $REFERENCE
wait

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ STAR1PASS-INDEX $ASSEMBLY TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "#@#############################################################"

echo
echo
echo

