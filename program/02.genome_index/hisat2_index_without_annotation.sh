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

hisat2-build -p $THREADS $REFERENCE \
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
