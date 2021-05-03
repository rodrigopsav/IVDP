#!/bin/bash

### Copy files to build snpEff database
mkdir -p $IVDP_DIR/dataBases/annotation_db/data/"$ASSEMBLY" > /dev/null 2>&1
cp $REFERENCE $IVDP_DIR/dataBases/annotation_db/data/"$ASSEMBLY"
cp $ANNOTATION $IVDP_DIR/dataBases/annotation_db/data/"$ASSEMBLY"

mv $IVDP_DIR/dataBases/annotation_db/data/"$ASSEMBLY"/*.fa $IVDP_DIR/dataBases/annotation_db/data/"$ASSEMBLY"/sequences.fa
mv $IVDP_DIR/dataBases/annotation_db/data/"$ASSEMBLY"/*.gtf $IVDP_DIR/dataBases/annotation_db/data/"$ASSEMBLY"/genes.gtf


### Add name of new database on "Non-Standard Databases" of snpEff.config file

LINE=$(cat $IVDP_DIR/dataBases/annotation_db/snpEff.config | grep -n -i "Non-Standard Databases" | awk -F":" '{print $1}')
LINE=$((LINE + 2))

awk -v line="$LINE" -v assembly="$ASSEMBLY" 'NR==line{print ""}1' $IVDP_DIR/dataBases/annotation_db/snpEff.config > $IVDP_DIR/dataBases/annotation_db/snpEff.config.tmp
mv $IVDP_DIR/dataBases/annotation_db/snpEff.config.tmp $IVDP_DIR/dataBases/annotation_db/snpEff.config
LINE=$((LINE + 1))
awk -v line="$LINE" -v assembly="$ASSEMBLY" 'NR==line{print "# "assembly}1' $IVDP_DIR/dataBases/annotation_db/snpEff.config > $IVDP_DIR/dataBases/annotation_db/snpEff.config.tmp
mv $IVDP_DIR/dataBases/annotation_db/snpEff.config.tmp $IVDP_DIR/dataBases/annotation_db/snpEff.config
LINE=$((LINE + 1))
awk -v line="$LINE" -v assembly="$ASSEMBLY" 'NR==line{print assembly".genome : "assembly}1' $IVDP_DIR/dataBases/annotation_db/snpEff.config > $IVDP_DIR/dataBases/annotation_db/snpEff.config.tmp
mv $IVDP_DIR/dataBases/annotation_db/snpEff.config.tmp $IVDP_DIR/dataBases/annotation_db/snpEff.config


### Build snpEff "Non-Standard Databases" database
snpEff -Xmx20g build -gtf22 -v $ASSEMBLY -c $IVDP_DIR/dataBases/annotation_db/snpEff.config
wait

