#!/bin/bash 

echo
echo "#@#############################################################"
echo "#@ DICTIONARY $ASSEMBLY"
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

picard CreateSequenceDictionary \
R=$REFERENCE \
O=$(echo "${REFERENCE%.fa}".dict)
wait

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ CREATE DICTIONARY $ASSEMBLY TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "#@#############################################################"

echo
echo
echo
