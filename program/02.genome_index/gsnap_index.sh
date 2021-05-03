#!/bin/bash

echo
echo "#@#############################################################"
echo "#@ GSNAP-INDEX $ASSEMBLY"
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

gmap_build -d $(echo $ASSEMBLY)idx -D $gsnap_index $REFERENCE
wait

cd $gsnap_index/"$(echo $ASSEMBLY)idx"/$(echo $ASSEMBLY)idx.maps
gtf_splicesites < $ANNOTATION | iit_store -o gsnap_splicesites.iit
gtf_introns < $ANNOTATION | iit_store -o gsnap_introns.iit

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ GSNAP-INDEX $ASSEMBLY TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo
echo "#@#############################################################"

echo
echo
echo

