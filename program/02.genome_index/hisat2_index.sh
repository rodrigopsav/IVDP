#!/bin/bash

echo
echo "#@#############################################################"
echo "#@ HISAT2-INDEX $ASSEMBLY"
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

extract_splice_sites.py -v \
$ANNOTATION \
> $hisat2_index/splice_sites.txt
wait

extract_exons.py -v \
$ANNOTATION \
> $hisat2_index/exons.txt
wait

hisat2-build -p $THREADS $REFERENCE \
--ss $hisat2_index/splice_sites.txt \
--exon $hisat2_index/exons.txt \
$hisat2_index/"$(echo $ASSEMBLY)idx"
wait


END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ HISAT2-INDEX $ASSEMBLY TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "#@#############################################################"

echo
echo
echo

