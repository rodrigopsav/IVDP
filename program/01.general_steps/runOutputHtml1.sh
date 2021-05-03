#!/bin/bash

# PAGES FOR HTML
# https://guides.codechewing.com/bash/reuse-html-template-bash
# https://www.w3schools.com/html/html_comments.asp
# https://html-color-codes.info/
# https://stackoverflow.com/questions/6913234/how-to-set-different-colors-in-html-in-one-statement

D1=$(date "+%D    %T")

######################################
### PARAMETERS
# Copy parameters
export IMAGE=$IVDP_DIR/program/01.general_steps/image
cp -r $IMAGE "$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/report/
cp $PARAMETERS "$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/report/"$(basename $PARAMETERS)"

sed -i 's/\\/\t/g' ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt
cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt | tail -n +2 >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/ivdp_sample_id.txt
sed -i $'1 i\\\nSAMPLE\tINDIVIDUAL_ID\tREAD_GROUP_ID' ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/ivdp_sample_id.txt
# Remove repeated header
cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/ivdp_sample_id.txt | uniq > ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/ivdp_sample_id.txt.bak
mv ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/ivdp_sample_id.txt.bak ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/ivdp_sample_id.txt
wait
# Remove repeated samples
cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/ivdp_sample_id.txt | head -n 1 > ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/ivdp_sample_id.txt.bak
cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/ivdp_sample_id.txt | tail -n +2 | sort | uniq >> ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/ivdp_sample_id.txt.bak
mv ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/ivdp_sample_id.txt.bak ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/ivdp_sample_id.txt
rm ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/"${ANALYSIS_ID}"_ivdp_sample_id.txt > /dev/null 2>&1

# Check errors in log folder
#"$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/errors.txt

# Remove directories
rm -r "$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/tmp_"${ANALYSIS_ID}" > /dev/null 2>&1
#rm -r "$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/chromosomes > /dev/null 2>&1
wait

if [[ $KEEP_LOG == "no" && ! -f "$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/errors.txt ]]; then
   rm -r "$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/log > /dev/null 2>&1
   rm -r "$OUTPUT_DIR"/ivdp_"$OUTPUT_NAME"/logSlurm > /dev/null 2>&1
fi
wait

# Count number of samples and individuals
export N_SAMPLES=$(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/ivdp_sample_id.txt | tail -n +2 | awk -F'\t' '{print $1}' | sort | uniq | wc -l)
export N_INDIVIDUALS=$(cat ${OUTPUT_DIR}/ivdp_${OUTPUT_NAME}/report/ivdp_sample_id.txt | tail -n +2 | awk -F'\t' '{print $2}' | sort | uniq | wc -l)

wait
######################################
# Simple automated HTML template

cat << _EOF_
<!DOCTYPE html>

<html>
   <style type="text/css">
   
   table, th, td {
   border: 1px solid black;
   }
   
   p.format1 { font-size:12px;line-height:50%;color:#C86904;font-weight:bold;font-family:verdana }
   span.format2 { font-size:12px;line-height:50%;color:#000000;font-weight:normal;font-family:verdana }
   pre.format3 { font-size:12px;line-height:50%;color:#086A87;font-weight:bold;text-align:center;font-family:verdana }
   p.format4 { font-size:12px;line-height:140%;color:#000000;font-weight:normal;text-align:center;font-family:verdana } 

   </style>

<head>
<TITLE> </TITLE>
</head>


<body style="background-color:#E0DCD1;margin-top:3px;margin-left:20px;margin-right:1px;" link="#0000FF" vlink="#0000FF" >

   <h2 style="font-size:30px;color:#336699;text-align:center;"><b>

<img src="report/image/DNA.png" class="center" width="40" height="40" />   </a>

					   IVDP - Integrated Variant Discovery Pipeline

   </b> </h2>


<!-- DETAILS OF THE ANALYSIS -->
   <p class="format1">DATE AND TIME: <span class="format2">$D1 </span> </p>
   <p class="format1">NAME OF ANALYSIS: <span class="format2">$OUTPUT_NAME </span> </p>
   <p class="format1">TYPE OF ANALYSIS: <span class="format2">$(if [[ $ANALYSIS == 1 ]]; then echo "WGS"; else echo "RNA-Seq"; fi) </span> </p>
   <p class="format1">TYPE OF FILES: <span class="format2">$(if [[ $TYPE == "se" ]]; then echo "Single-end"; else echo "Paired-end"; fi) </span> </p>
   <p class="format1">NAME OF REFERENCE GENOME: <span class="format2">$ASSEMBLY </span> </p>
   <p class="format1">GENE COUNT PROGRAMS: <span class="format2">$GENECOUNT_PROGRAM </span> </p>
   <p class="format1">ALIGNMENT PROGRAMS: <span class="format2">$ALIGNER_PROGRAM </span> </p>
   <p class="format1">VARIANT CALLING PROGRAMS: <span class="format2">$CALLER_PROGRAM </span> </p>
   <p class="format1">MINIMUM READ LENGTH: <span class="format2">$MIN_READ_LENGTH </span> </p>
   <p class="format1">MAXIMUM READ LENGTH: <span class="format2">$MAX_READ_LENGTH </span> </p>
   <p class="format1">MINIMUN DEPTH FOR FILTERED LOCI: <span class="format2">$MIN_DEPTH </span> </p>
   <p class="format1">MINIMUN ALLELE FREQUENCY FOR FILTERED LOCI: <span class="format2">$MAF </span> </p>
   <p class="format1">MAXIMUN GENOTYPES MISSING RATE: <span class="format2">$MISSING </span> </p>
<!--   <p class="format1">CHROMOSOME INFORMATION: <span class="format2"><a href="chromosomes">click here </a> </span> </p>  -->
   <p class="format1">NUMBER OF FASTQ FILES (SAMPLES): <span class="format2">$N_SAMPLES </span> </p>
   <p class="format1">NUMBER OF INDIVIDUALS: <span class="format2">$N_INDIVIDUALS </span> </p>
   <p class="format1">LIST OF SAMPLES: <span class="format2"><a href="report/ivdp_sample_id.txt">click here </a> </span> </p>
   <p class="format1">PARAMETER FILE: <span class="format2"><a href="report/$(basename $PARAMETERS)">click here </a> </span> </p>

<br />

<table>
<!-- TABLE HEADER -->
  <tr>
    <th>DATA FILES</th>
    <th>REPORTS</th>
    <th>STATISTICS</th>
  </tr>

<!-- FIRST COLUMN -->     
   <tr> <td rowspan="15">				      						      
     <pre class="format3">  QUALITY CONTROL OF FASTQ FILES  </pre>
     <p class="format4">
<a href="data/fastqQC">QC AND TRIMMED FASTQ FILES</a><br />
<br />
     </p>
     
     <pre class="format3">  ALIGNMENT AND GENE COUNTS  </pre>
     <p class="format4"> 
<a href="data/alignment">ALIGNED BAM FILES</a><br />
<a href="data/gene_count">GENE COUNTS</a><br />
<br />
     </p>
   
     <pre class="format3">  VARIANTS  </pre>
     <p class="format4">    
<a href="data/gvcf">GVCF FILES</a><br />
<a href="data/vcf/raw">RAW VCF FILES</a><br />
<a href="data/vcf/filtered">FILTERED VCF FILES</a><br />
<a href="data/vcf/combined">COMBINED SNP VCF FILES</a><br />
     </p>
   </td> </tr>



<!-- SECOND COLUMN -->
   <tr> <td rowspan="15">				      						      
     <pre class="format3">  QUALITY CONTROL OF FASTQ FILES  </pre>
     <p class="format4">
<a href="report/fastqQC"> QC AND TRIMMED FASTQ FILES </a><br />
<br />
     </p>
     
     <pre class="format3">  ALIGNMENT AND GENE COUNTS  </pre>
     <p class="format4"> 
<a href="report/alignment">ALIGNED BAM FILES</a><br />
<a href="report/gene_count">GENE COUNTS</a><br />
<br />
     </p>
     
     <pre class="format3">  VARIANTS  </pre>
     <p class="format4">
<br />
<a href="report/vcf/raw">RAW VCF FILES</a><br />
<a href="report/vcf/filtered">FILTERED VCF FILES</a><br />
<a href="report/vcf/combined">COMBINED SNP VCF FILES</a><br />
     </p>
   </td> </tr>     



<!-- THIRD COLUMN --> 
   <tr> <td rowspan="15">
     <pre class="format3">  QUALITY CONTROL OF FASTQ FILES  </pre>
     <p class="format4">
<a href="stats/fastqQC">QC AND TRIMMED FASTQ FILES</a><br />
<br />
     </p>

     <pre class="format3">  ALIGNMENT  </pre>
     <p class="format4">    
<a href="stats/alignment">ALIGNMENT</a><br />
<br />
<br />
     </p>
   
     <pre class="format3">  VARIANTS  </pre>
     <p class="format4">
<br />
<a href="stats/vcf/raw">RAW VCF FILES</a><br />
<a href="stats/vcf/filtered">FILTERED VCF FILES</a><br />
<a href="stats/vcf/combined">COMBINED SNP VCF FILES</a><br />
     </p>
   </td> </tr>
   
</table>   
<p><strong>Note:</strong> Some directories can be empty depending on the options used in the parameter file.</p>
   
<address>
<br />
Written by <a href="mailto:rodrigopsa@yahoo.com">Rodrigo Savegnago</a>.<br> 
<a href="https://github.com/rodrigopsav/IVDP">  <img src="report/image/Github.png" width="35" height="35" />   </a>
<a href="https://www.linkedin.com/in/rodrigopsav/">  <img src="report/image/Linkedin.png" width="35" height="35" />   </a>
<a href="https://www.facebook.com/rodrigopsav/"> <img src="report/image/Facebook.png" width="35" height="35" />   </a>
<a href="https://www.instagram.com/rodrigopsav/">  <img src="report/image/Instagram.png" width="35" height="35" />   </a>
<a href="https://twitter.com/rodrigopsav">  <img src="report/image/Twitter.png" width="35" height="35" />   </a><br>

Michigan State University<br>
Department of Animal Science<br>
474 S Shaw Ln, East Lansing, MI<br>
48824, USA
</address>

</body>

</html>
_EOF_

