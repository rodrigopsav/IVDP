
# IVDP: Integrated Variant Discovery Pipeline
## User manual and guide

IVDP is a collection of Bash and R scripts developed for calling variants purposes - SNPs (single-nucleotide polymorphisms) and Indels (insertions and deletions) - from Whole Genome Sequence (WGS) and RNA Sequence (RNAseq) data. It can also do gene counts of RNAseq data.
IVDP combines more than 30 programs and R packages for data manipulation and bioinformatic analysis. It is optimized to run in a local server machine or HPCC with slurm scheduler.


[Citation](#citation)   
[Running IVDP first time: quick steps](#running-ivdp-first-time-quick-steps)   
[Download IVDP](#download-ivdp)   
[Install IVDP dependencies](#install-ivdp-dependencies)   
[Using conda environments without IVDP](#using-conda-environments-without-ivdp)   
[Running IVDP](#running-ivdp)   
[Killing IVDP analysis](#killing-ivdp-analysis)   
[Checking main log](#checking-main-log)   
[IVDP parameter file](#ivdp-parameter-file)   
[Slurm parameter file](#slurm-parameter-file)   
[IVDP Output](#ivdp-output)   
[IVDP Examples](#ivdp-examples)   
* [Example 1](#example-1): WGS analysis   
* [Example 2](#example-2): WGS analysis   
* [Example 3](#example-3): WGS analysis   
* [Example 4](#example-4): WGS analysis   
* [Example 5](#example-5): WGS analysis   
* [Example 6](#example-6): WGS analysis   
* [Example 7](#example-7): Paired-end RNAseq analysis   
* [Example 8](#example-8): Paired-end RNAseq analysis   
* [Example 9](#example-9): Paired-end RNAseq analysis   
* [Example 10](#example-10): Single-end RNAseq analysis   
* [Example 11](#example-11): Single-end RNAseq analysis   
* [Example 12](#example-12): Single-end RNAseq analysis   
* [Example 13](#example-13): Single-end RNAseq analysis   

[Caveats](#caveats)   
[Troubleshooting](#troubleshooting)   


## Citation
Savegnago, R. P. IVDP: Integrated Variant Discovery Pipeline for whole genome sequencing and RNA sequencing. 2021. GitHub repository, https://github.com/rodrigopsav/IVDP

[Download citation](https://github.com/rodrigopsav/IVDP/tree/main/program/01.general_steps/image/ivdp_citation.bib)

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

## Running IVDP first time: quick steps

```
# Download IVDP
git clone https://github.com/rodrigopsav/IVDP.git
wait

# Install IVDP dependencies
cd IVDP/install_ivdp_dependencies
./install_ivdp_dependencies.sh -d .
cd ..
wait

# Running IVDP in a local server (pe_ex1a.txt)
./ivdp.sh -p params/pe_ex1a.txt

# Running IVDP in HPCC with slurm (pe_ex1a.txt)
# ./ivdp.sh -p params/pe_ex1a.txt -c params/configSlurm.txt
```

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

## Download IVDP

```
git clone https://github.com/rodrigopsav/IVDP.git 
```

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

## Install IVDP dependencies

IVDP dependencies are installed in three conda environments (ivdp, ivdp2, and r-env4.0). Before run IVDP, you **MUST** install the dependencies (IVDP/install_ivdp_dependencies/install_ivdp_dependencies.sh file) even if you have already installed the programs in your machine. To install IVDP dependencies, run:
```
cd IVDP/install_ivdp_dependencies
./install_ivdp_dependencies.sh -d <directory/to/install/ivdp/dependencies>
```
The bash script (./install_ivdp_dependencies.sh) will install conda (if it has not installed before) and the programs required by IVDP. 

Unfortunately, incompatibilities can happen when install the most recent version of the programs. If you detect some errors after install IVDP dependencies, use the following command lines to remove the previous installation and re-install them with the versions used originally to develop IVDP:
```
conda env remove --name ivdp
conda env remove --name ivdp2
./install_ivdp_dependencies_versions.sh -d <directory/to/install/ivdp/dependencies>
```

Please visit the official websites of the programs used on IVDP to learn more about them:

--**Bioinformatic programs**: [BBMap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/)	  [BCFtools](http://www.htslib.org/download/)	  [BEDtools](https://bedtools.readthedocs.io/en/latest/)	  [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml)   [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)	  [BWA-MEM](https://github.com/lh3/bwa)	  [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2)    [GNU-datamash](https://www.gnu.org/software/datamash/)    [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)    [featureCounts](http://subread.sourceforge.net/)    [Freebayes](https://github.com/ekg/freebayes)   [Freebayes-parallel](https://bioconda.github.io/recipes/freebayes/README.html)    [GATK](https://github.com/broadinstitute/gatk)    [GSNAP](http://research-pub.gene.com/gmap/)   [HISAT2](https://daehwankimlab.github.io/hisat2/)   [HTSeq](https://htseq.readthedocs.io/en/master/)    [HTSlib](http://www.htslib.org/download/)   [MultiQC](https://multiqc.info)   [Picard](https://broadinstitute.github.io/picard/)    [Platypus](https://github.com/andyrimmer/Platypus)    [Plink1.9](https://www.cog-genomics.org/plink/)    [Plink2](https://www.cog-genomics.org/plink/2.0/)    [Sambamba](https://lomereiter.github.io/sambamba/index.html)    [Samtools](http://www.htslib.org/download/)   [snpEff and snpSift](http://snpeff.sourceforge.net)   [STAR](https://github.com/alexdobin/STAR)   [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)    [vcfstats](https://pypi.org/project/vcfstats/)   [VCFtools](https://vcftools.github.io/index.html)   [Zlib](http://www.zlib.net)

--**Programs for file manipulation**: [dos2unix](https://pypi.org/project/dos2unix/)    [ghostscript](https://www.ghostscript.com/)    [img2pdf](https://pypi.org/project/img2pdf/)    [pigz](http://zlib.net/pigz/)    [sed](https://pypi.org/project/sed/)    [vcflib](https://github.com/vcflib/vcflib#INSTALL)    [zlib](http://zlib.net/)

-- **R packages**: [bigmemory.sri](https://cran.r-project.org/web/packages/bigmemory.sri/index.html)    [CMplot](https://github.com/YinLiLin/CMplot)    [data.table](https://cran.r-project.org/web/packages/data.table/index.html)    [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)    [formattable](https://cran.r-project.org/web/packages/formattable/index.html)    [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)    [gridextra](https://cran.r-project.org/web/packages/gridExtra/index.html)        [htmltools](https://cran.r-project.org/web/packages/htmltools/index.html)    [rmvp](https://github.com/xiaolei-lab/rMVP)    [webshot](https://cran.r-project.org/web/packages/webshot/index.html)

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

## Using conda environments without IVDP

You can use all the programs without IVDP, by activating the proper conda environment:

```bash
conda activate ivdp
conda activate ivdp2
conda activate r-env4.0
```
Only freebayes and platypus are installed in ivdp2 environment with python 2.7 while the others are installed in ivdp conda environment with python 3.8. The R program v4 and its packages are in r-en4.0. To deactivate a conda environment, use:
```bash
conda deactivate
```
It is usual to activate only one conda environment at time to avoid incompatibilities. However, if you want to use more than one at once, type:
```bash
conda activate ivdp
conda activate --stack ivdp2
```
**NOTE**: If you activate two conda envs at the same time, some programs cannot work properly.

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

## Running IVDP

**Run ./ivdp.sh from IVDP folder to avoid possible problems.**

To run IVDP in a local machine, type:
```
./ivdp.sh -p <ivdp parameter file>
```

and to run on a HPCC with slurm scheduler, type:
```
./ivdp.sh -p <ivdp parameter file> -c <slurm config file>
```

The IVDP parameter file is mandatory.  The slurm config file will activate the jobs submission on HPCC with slurm.

**NOTE 1**: IVDP parameter file and slurm config file **MUST BE** anywhere in IVDP folder.

**NOTE 2**: Take care to set up the directories path in ivdp parameter file **before** running IVDP.

**NOTE 3**: Do not use nohup and & to run IVDP in the background. After ./ivdp.sh starts, all the processes will be automatically sent to the backgroud (in a local machine).

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

## Killing IVDP analysis
In the end of each submission, IVDP shows a message to kill the analysis like that:

```
# IVDP running in a local server
To kill this IVDP analysis, run: kill -- 1463

# IVDP running on HPCC with slurm
for job in $(squeue -u $USER | grep 463739 | awk '{print \$1}'); do scancel $job; done
```

**NOTE**: each analysis shows a different number (like 1463 and 463739 in the example above. Copy and paste these command lines anywhere until  the analysis finishes (just in case you want to stop it earlier).

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

## Checking main log
To check the main log for **non HPCC analysis**, type:
```
./ivdpLog.sh -l 1 -a <analysis ID> -o <ivdp output dir>
```
To skip it anytime, press Ctrl+C. 

To access the IVDP log for t**HPCC** jobs, type:
```
./ivdpLog.sh -l 2 -a <analysis ID> -o <ivdp output dir>
```
**NOTE**: IVDP shows these command lines with \<analysis ID> and \<output directory> in the end of each submission.

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

## IVDP parameter file
IVDP parameter file contains variables used by IVDP:

**########## Input data section ##########**

* **INPUT_DIR**: Directory path with the original fastq files.
  > **NOTE**: Do not use ~/ to indicate home directory. Use the absolute path: /**home**/USERNAME/

* **OUTPUT_DIR**: Output directory folder (Default is the home directory if left blank: OUTPUT_DIR=).
  > **NOTE**: Do not use ~/ to indicate home directory. Use the absolute path: /**home**/USERNAME/

* **OUTPUT_NAME**: Analysis name. Choose a single name without spaces (Default is user name if left blank: OUTPUT_NAME=).

* **LIST_SAMPLES**: Text file path with a list of samples id and individual id, separated by space [(go to Example 3 for more details)](#example-3). Use LIST_SAMPLES=none if no list is available.
  > **NOTE**: Do not use ~/ to indicate home directory. Use the absolute path: /**home**/USERNAME/

* **ANALYSIS**: Choose between Whole Genome Sequence (**ANALYSIS=1**) or RNAseq data (**ANALYSIS=2**).

* **EXTENSION_SE**: For single-end fastq files only. Write the extension of original fastq files. If paired-end files, left blank (EXTENSION_SE=) or remove from parameter file.

* **EXTENSION_PE1**: For paired-end fastq files only. Write the extension of PE1 fastq files. If single-end files, left blank (EXTENSION_PE1=) or remove from parameter file.

* **EXTENSION_PE2**: For paired-end fastq filesonly. Write the extension of PE2 fastq files. If single-end files, left blank (EXTENSION_PE2=) or remove from parameter file.

**########## Input genome section ##########**

* **ASSEMBLY**: Assign a name for the reference genome. Eg. ASSEMBLY=hg38 (**Do not use spaces!**).

* **REFERENCE**: File directory for reference genome **(*.fa)**. The reference genome must be **uncompressed file** and not **(\*.fa.gz)**.
  > **NOTE**: Do not use ~/ to indicate home directory. Use the absolute path: /**home**/USERNAME/

* **VARIANTS**: File directory for variants file **(\*.vcf.gz)** for the reference genome. Use VARIANTS=none if no variants is available.The variants must be **(\*.vcf.gz)** and not **(\*.vcf)**.
  > **NOTE**: Do not use ~/ to indicate home directory. Use the absolute path: /**home**/USERNAME/

* **ANNOTATION**: File directory for annotation file **(\*.gtf)** for the reference genome. Use ANNOTATION=none if no annotation is available. The annotation must be **uncompressed file** and not **(\*.gtf.gz)**.
  > **NOTE**: Do not use ~/ to indicate home directory. Use the absolute path: /**home**/USERNAME/

**########## Workflow steps section ##########**

* **STEP_QC_TRIM** Activate quality control and adapter trimming step? STEP_QC_TRIM=yes or STEP_QC_TRIM=no (Default yes).

* **STEP_ALIGNMENT**: Activate alignment step? STEP_ALIGNMENT=yes or STEP_ALIGNMENT=no (Default yes).

* **STEP_GENECOUNT**: # Activate gene count step? STEP_GENECOUNT=yes or STEP_GENECOUNT=no (Default: yes for RNAseq and no for WGS).

* **STEP_MDUP**: Activate mark duplicated read step? STEP_MDUP=yes or STEP_MDUP=no (Default: yes). It includes splitncigar step for RNAseq.

* **STEP_BQSR**: Activate base quality score recalibration step? STEP_BQSR=yes or STEP_BQSR=no (Default: yes). Only works if a valid anotation.vcf file is supplied (**VARIANTS** variable of IVDP parameter file)

* **STEP_VAR_CALL**: Activate variant calling step? STEP_VCF_CALL=yes or STEP_VCF_CALL=no (Default: yes).

* **STEP_GVCF_TO_VCF**: Convert gvcf to vcf? STEP_GVCF_TO_VCF=yes or STEP_GVCF_TO_VCF=no (Default: yes). **NOTE**: This option works only with specific variant callers (gatk4GenomicsDBImport or gatk4CombineGVCFs).

* **STEP_VCF_FILTER**: Activate vcf filtering step? STEP_VCF_FILTER=yes or STEP_VCF_FILTER=no (Default: yes).

**########## Choose the programs ##########**

* **ALIGNER_PROGRAM**: Choose programs for alignment step. If more than one program are chosen, use ',' between them (Options: bbmap, bowtie, bowtie2, bwa, bwa2, gsnap, hisat2, star).
  > **NOTE1**: If ALIGNER_PROGRAM is empty (ALIGNER_PROGRAM=) or if it is not in the parameter file, the default is **ALIGNER_PROGRAM=bwa2** if ANALYSIS=1 and **ALIGNER_PROGRAM=hisat2** if ANALYSIS=2.

* **CALLER_PROGRAM**: Choose programs for variant calling step. If more than one program are chosen, use ',' between them (Options: bcftools, freebayes, freebayes-parallel, gatk4, gatk4GenomicsDBImport, gatk4CombineGVCFs, platypus, varscan).
  > **NOTE1**: If CALLER_PROGRAM is empty (CALLER_PROGRAM=) or if it is not in the parameter file, the default is **CALLER_PROGRAM=gatk4GenomicsDBImport** if ANALYSIS=1 and **CALLER_PROGRAM=gatk4** if ANALYSIS=2.
  
  > **NOTE2**: freebayes and freebayes-parallel are the same program, but the second one is faster.
  
  > **NOTE3**: gatk4 calls the variants using GATK HaplotypeCaller from a multi-sample the [joint calling approach](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants). The gatk4GenomicsDBImport and gatk4CombineGVCFs options use GATK HaplotypeCaller with gVCF model, i.e. a gVCF file will be created for each individual and after that merged into a single file using [GATK GenomicsDBImport](https://gatk.broadinstitute.org/hc/en-us/articles/360057439331-GenomicsDBImport) or [GATK CombineGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360056969772-CombineGVCFs). Both options do the same job. Do not use both at the same time. Generally, gatk4GenomicsDBImport is faster than gatk4CombineGVCFs.

   > **NOTE3**: Use gatk4CombineGVCFs if gatk4GenomicsDBImport fails or if you have repeated samples from the same individual (i.e. various fastq files for same individual). Go to [Example 3](#example-3) to understand how to handle with repeated runs.

* **GENECOUNT_PROGRAM**: Choose programs for gene count step. If more than one program are chosen, use ',' between them (Options: htseq, featurecounts). Use GENECOUNT_PROGRAM=none if you do not want gene counts.
  > **NOTE1**: GENECOUNT_PROGRAM only works with RNAseq analysis (**ANALYSIS=2**). This option will be deactivated with WGS (**ANALYSIS=1**) even if STEP_GENECOUNT=yes and GENECOUNT_PROGRAM is different of "none".
  
  > **NOTE2**: GENECOUNT_PROGRAM requires a valid annotation gft file (**ANNOTATION** variable of IVDP parameter file)

* **FEATURE_TYPE**: Any feature of the third column of gtf file. (Default FEATURE_TYPE=exon). Use FEATURE_TYPE=none if you do not want gene counts. To check the possible features on third column of the gtf file, type on linux terminal:

   ```bash
    $ grep -v "#" file.gtf | awk '{print $3}' | sort | uniq
    ```

  > **NOTE1**: FEATURE_TYPE only works with RNAseq analysis (**ANALYSIS=2**). This option will be deactivated with WGS (**ANALYSIS=1**) even if STEP_GENECOUNT=yes.

**########## General parameters section ##########**

* **CALL_BY_CHROM**: Should variant calling be splitted by chromosome? CALL_BY_CHROM=no or CALL_BY_CHROM=yes (Default yes).
  > **NOTE1**: it works only on local servers.
 
* **BP_BY_CHROM**: This option breaks the chromosomes in smaller chunks to speed up the variant calling process (e.g. BP_BY_CHROM=20000 breaks each chromosome in chunks of 20mb). BP_BY_CHROM=all does not break the chromosome (default). 
  > **NOTE1**: it works only on HPCC.

  > **NOTE2**: If the variable is empty (BP_BY_CHROM=) or if it is removed from parameter file, IVDP will use the default.
 
 * **CHROM_SET**: Number of autosome pairs of the specie.

* **SELECT_SET**: select which chromosomes you wish in vcf files. (e.g. SELECT_SET=1,2,3). SELECT_SET=all include all autosomes (default). If variable is empty (SELECT_SET=) or if it is removed from parameter file, IVDP will use the default.

* **MIN_READ_LENGTH**: Minimum read length to be kept after adapter trimming (Default 36).

* **MAX_READ_LENGTH**: Maximum read length to be kept after adapter trimming (Default 150).

* **DOWNSAMPLING_FASTQ**: Down sampling the fastq files according to the depth of coverage chosen. (Default 0). DOWNSAMPLING_FASTQ=0 means no reads sampling from fastq files.  DOWNSAMPLING_FASTQ=1 means downsampling reads from fastq files to 1x of depth of coverage.

   > **NOTE:** Down-sampling reads can save computational resources and analysis time when the fastq files are extremely large. For example, it is possible to down-sampling the reads to decrease the coverage from 30x to 15x to call the variants. The relation between the number of reads, average read length and coverage can be calculated using:

   C = (N\*L)/G, which C=coverage, N=number of reads, L=average read length (in base pairs), and G=size of reference genome (in base pairs). **Multiply this formula by 2 if fastq files are paired-end**.

   > Considering a paired-end fastq file (multiply the formula by 2) with 450 million of reads in each file, average read length of 100bp, and the human reference genome (which contains approximately 3 billion of base pairs), the coverage will be approximately: 
   C = (450,000,000 \* 100 \* 2) / 3,000,000,000 = 30x

   > Downsampling the fastq file from 30x to 15x of coverage:
   15 = (N \* 100 \* 2) /3,000,000,000 --> N = 225 million reads.
   
   > In this case, if you choose  **DOWNSAMPLING_FASTQ=15** IVDP will random sampling 225 million reads from the fastq files.

* **DOWNSAMPLING_BAM**: The same idea as **DOWNSAMPLING_FASTQ**, but downsampling reads from the raw bam file.
  > **NOTE:**Do not use DOWNSAMPLING_FASTQ and DOWNSAMPLING_BAM at the same time. Choose only one.
 
* **MIN_DEPTH**: Minimum locus depth of coverage for the filtering step (Default 3, if STEP_VCF_FILTER=yes).

* **MAF**: Minimum allele frequency for the filtering step (Default 0.005, if STEP_VCF_FILTER=yes).

* **MISSING**: Maximum genotype-based missing rate to keep a locus (Default 0.3, if STEP_VCF_FILTER=yes). MISSING=0 means no missing values are allowed while MISSING=1 means no filter for missing rate. 

* **COMBINE_VCF**: Choose COMBINE_VCF=none, COMBINE_VCF=partial, COMBINE_VCF=full, or COMBINE_VCF=partial,full. For more details, go to [Example 5](#example-5)

  > ***COMBINE_VCF=none:*** do not combine different vcfs;
  > 
  > ***COMBINE_VCF=partial:*** combine common SNPs that appeared AT LEAST in two different vcf files if more than one variant caller is used;
  > 
  > ***COMBINE_VCF=full:*** combine common SNPs that appeared in ALL the vcf files if if more than one variant caller is used.

  > ***NOTE1:*** How COMBINE_VCF works?
  > Let's suppose you chose three variant callers (i.e. bcftools, freebayes and gatk4). In the end, IVDP will output one VCF file from each variant caller. If you use **COMBINE_VCF=partial**, IVDP will combine these three VCFs in a single one by selecting the SNPs that appeared at least on two different VCF files. If you choose **COMBINE_VCF=full**, IVDP will combine these three VCFs in a single one selecting only the SNPs that appeared in all the VCF files. 
  > Now let's suppose you chose two aligners (bowtie2 and bwa) and three variant callers (bcftools, freebayes and gatk4). In the end, IVDP will output six VCF files (bowtie2_bcftools, bowtie2_freebayes, bowtie2_gatk4, bwa_bcftools, bwa_freebayes, and bwa_gatk4). No matter if you chose **COMBINE_VCF=partial** or **COMBINE_VCF=full**, IVDP will combine the VCFs within aliners. In this case, IVDP will output two VCF files: one VCF file from the combination among bcftools, freebayes and gatk4 from **bowtie2** and another VCF file from the combination among bcftools, freebayes and gatk4 from **bwa**. 

  > **COMBINE_VCF** works only if more than one variant caller is used **AND** if STEP_VCF_FILTER=yes. IVDP combined the filtered vcfs and not the raw vcfs.

* **THREADS**: Number of threads for analysis.
  > **NOTE1**: it works only on local servers.

* **BATCH**: Number of samples to be processed at same time.
  > **NOTE1**: it works only on local servers.

* **STATS_ALIGNMENT**: Should IVDP calculate stats from alignment step? STATS_ALIGNMENT=no or STATS_ALIGNMENT=yes (Default yes)

* **KEEP_LOG**: Keep the log files in the output folder? KEEP_LOG=no or KEEP_LOG=yes (Default yes)

* **KEEP_INTERMEDIATE_DATA**: Keep intermediate fastq and bam files in the output folder? KEEP_INTERMEDIATE_DATA=no or KEEP_INTERMEDIATE_DATA=yes (Default yes)

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

## Slurm parameter file
A  text file with slurm parameters  is required to activate the jobs submissions on HPCC with slurm scheduler. There are 3 variables **for each** IVDP step:
1) time required for the analysis;
2) number of cores;
3) amount of ram memory.

For example, to trimming the raw fastq files, I chose 8 hours to complete the jobs, a cpu with 16 cores and 32G of ram memory:

TRIMMING_TIME=08:00:00   
TRIMMING_CPU=16   
TRIMMING_MEM=32G   


All the variables are in slurm parameter file are:

**Fastq QC and Trimming**   
TRIMMING_TIME=08:00:00   
TRIMMING_CPU=16   
TRIMMING_MEM=32G   
   
**Alignment**   
ALIGNMENT_TIME=48:00:00   
ALIGNMENT_CPU=16   
ALIGNMENT_MEM=40G   
   
**Gene counts**   
GENECOUNT_TIME=04:00:00   
GENECOUNT_CPU=16   
GENECOUNT_MEM=32G   
   
**Mark duplicated reads**   
MDUP_TIME=04:00:00   
MDUP_CPU=8   
MDUP_MEM=32G   
   
**Base quality score recalibration**   
BQSR_TIME=04:00:00   
BQSR_CPU=8   
BQSR_MEM=32G   
   
**Variant calling**   
CALLVARIANT_TIME=96:00:00   
CALLVARIANT_CPU=16   
CALLVARIANT_MEM=40G   
   
**Convert gvcf to vcf**   
GVCF_TO_VCF_TIME=96:00:00   
GVCF_TO_VCF_CPU=8   
GVCF_TO_VCF_MEM=64G   
   
**VCF filtering**   
FILTERVCF_TIME=04:00:00   
FILTERVCF_CPU=4   
FILTERVCF_MEM=16G   
   
The amount of resource can vary depending on the programs used and the amout of data and the quality of the  data (depth of coverage, genotype missing rate, etc). Each step will be activated or not depending on the options in the main ivdp parameter file. So, you can specify all these slurm parameters, but they will be used depending on the options in the main parameter file.

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>


## IVDP Output
IVDP outputs a series of files, statistics and reports organized in four folders: data, log, report, and stats. Click on the html file to access these folders and subfolders in a more intuitive way. The html output file looks like:

![alt text](https://github.com/rodrigopsav/IVDP/blob/main/program/01.general_steps/image/html_output.png)

The first lines in orange show a brief description of the analysis. You can also access a list with samples and the ivdp parameter file.
The table contain three columns:
<br>

#### FIRST COLUMN: DATA FILES
 
 * <ins>QUALITY CONTROL OF FASTQ FILES</ins>: it contains the fastq.gz files after adapter removal and quality control by trimmomatic.
    
  * <ins>ALIGNMENT AND GENE COUNTS</ins>: it contains raw bam files (with read group and sorted by position), mdup (with mark duplicated reads), and bqsr (mdup + base quality score recalibration). The gene counts will be available only for RNAseq data (***ANALYSIS***=2 in ivdp parameter file).
    
  * <ins>VARIANTS</ins>:  It contains the gvcfs (only for ***CALLER_PROGRAM***=gatk4GenomicsDBImport or ***CALLER_PROGRAM***=gatk4CombineGVCFs), the raw vcf with rsID (if VARIANTS is supplied in ivdp parameter file) and annotated with [snpEff and snpSift](http://snpeff.sourceforge.net) (if ***ANNOTATION*** is supplied in ivdp parameter file) , filtered vcf (if ***STEP_VCF_FILTER***=yes), and combined vcf (if more than one variant caller is chosen: eg. ***CALLER_PROGRAM***=bcftools, gatk4).
  
    > **NOTE1**: IVDP will remove trimmed fastq, raw and mdup bam files if ***KEEP_INTERMEDIATE_DATA***=no

     > **NOTE2**: if a first analysis starts with ***STEP_QC_TRIM***=no, the trimmed fastq files will not be created and IVDP will use the raw fastq files for the alignment step. <br>

#### SECOND COLUMN: REPORTS
It contains html and PDF reports from files in the first column.

> **NOTE1**: the *zip file in the IVDP output folder contain all  reports that can be saved in other directory.

> **NOTE2**: the alignment report contain 3 informations among others: depth of coverage, breadth of coverage and depth of coverage only for covered regions.
A brief explanation about depth and breadth of  coverage can be found [here](https://sites.google.com/site/wiki4metagenomics/pdf/definition/coverage-read-depth). The depth of coverage only for covered regions is a conditional coverage: given a percentage of the genome that was covered genome by reads, what is the depth of coverage? This is the concept of depth of coverage only for covered regions.

#### THIRD COLUMN: STATISTICS
  
  * <ins>QUALITY CONTROL OF FASTQ FILES</ins>
    - *fastq_raw:* fastqc statistics for each original sample.
    - *fastq_trimmed:* fastqc statistics for each trimmed sample.
    - *trimmomatics:* statistics from trimmomatic.

  * <ins> ALIGNMENT </ins>: 
     * .coverageBed: depth coverage for each locus from [BEDtools](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html).
     * .propcoverage: proportion of coverage from [BEDtools](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html)

  * <ins> VARIANTS </ins>: for raw, filter and combined
    * .afreq: allele frequencies from [plink](https://www.cog-genomics.org/plink/2.0/formats#afreq)
    * .het: inbreeding coefficient from [plink](https://www.cog-genomics.org/plink/2.0/formats#het)
    * .smiss: sample-based missing rate from [link](https://www.cog-genomics.org/plink/2.0/formats#smiss)
    * .vmiss: genotype-based missing rate from [plink](https://www.cog-genomics.org/plink/2.0/formats#vmiss)
    * .stats: vcf statistics from [bcftools stats](http://samtools.github.io/bcftools/bcftools.html#stats)
    * snpEff_effects and snpEff_counts: statistics from [snpEff](https://pcingola.github.io/SnpEff/se_outputsummary/)

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>


## IVDP Examples

### Whole Genome Sequencing (examples 1 to 5)

**NOTE**: Do not forget to change the directories in IVDP parameter files before running the examples.

#### Example 1

**a. Standard WGS analysis**: WGS (ANALYSIS=1) with one aligner (bwa2), and one variant caller (gatk4GenomicsDBImport). If IVDP is in a local server, 5 samples will run at the same time (BATCH=5) with 40 cores (THREADS=40). It it runs on HPCC, all the samples will be submitted in parallel at once.
```bash
# Run IVDP on a local server
./ivdp.sh -p params/pe_ex1a.txt
# Run IVDP on HPCC with slurm
./ivdp.sh -p params/pe_ex1a.txt -c params/configSlurm.txt
``` 

**b.**:  Same as 1a, but **downsampling fastq** files to 10x of depth of coverage (**DOWNSAMPLING_FASTQ=10**).
```bash
# Run IVDP on a local server
./ivdp.sh -p params/pe_ex1b.txt
# Run IVDP on HPCC with slurm
./ivdp.sh -p params/pe_ex1b.txt -c params/configSlurm.txt
```

**c.**: Same as 1a, but **downsampling bam** files to 10x of depth of coverage (**DOWNSAMPLING_BAM=10**).
```bash
# Run IVDP on a local server
./ivdp.sh -p params/pe_ex1c.txt
# Run IVDP on HPCC with slurm
./ivdp.sh -p params/pe_ex1c.txt -c params/configSlurm.txt
```

**d.**: Same as 1b, but skipping read trimming step (**STEP_QC_TRIM=no**). This will align the raw fastq without applying any quality control. It is useful when the fastq were filtered previously.
```bash
# Run IVDP on a local server
./ivdp.sh -p params/pe_ex1d.txt
# Run IVDP on HPCC with slurm
./ivdp.sh -p params/pe_ex1d.txt -c params/configSlurm.txt
```

**e.**: Same as 1c, but skipping read trimming step (**STEP_QC_TRIM=no**).
```bash
# Run IVDP on a local server
./ivdp.sh -p params/pe_ex1e.txt
# Run IVDP on HPCC with slurm
./ivdp.sh -p params/pe_ex1e.txt -c params/configSlurm.txt
```

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

#### Example 2
Same as example 1a, but using bowtie2 as aligner and **two variant callers** (bcftools and gatk4GenomicsDBImport). The 2 filtered vcfs will be combined into a single vcf (**COMBINE_VCF=partial**) keeping the common positions, reference and alternate alleles.
```bash
# Run IVDP on a local server
./ivdp.sh -p params/pe_ex2.txt
# Run IVDP on HPCC with slurm
./ivdp.sh -p params/pe_ex2.txt -c params/configSlurm.txt
```

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

#### Example 3
**a.**:  Same as example 1a, but using a list of samples to subset 3 individuals out of 5 (**LIST_SAMPLES=/home/work/rps/softwares/IVDP/examples/list1.txt**) in folder examples/pe_wgs2.
```bash
# Run IVDP on a local server
./ivdp.sh -p params/pe_ex3a.txt
# Run IVDP on HPCC with slurm
./ivdp.sh -p params/pe_ex3a.txt -c params/configSlurm.txt
```

**b.**:  Same as example 3a, but considering repeated runs for the same individual (**LIST_SAMPLES=/home/work/rps/softwares/IVDP/examples/list2.txt**) in folder examples/pe_wgs2.

> **NOTE1**: Always use gatk4CombineGVCFs when call variants from repeated runs of the same sample!

```bash
# Run IVDP on a local server
./ivdp.sh -p params/pe_ex3b.txt
# Run IVDP on HPCC with slurm
./ivdp.sh -p params/pe_ex3b.txt -c params/configSlurm.txt
```

> **NOTE2: How to create a list of samples?**
The list of samples consists of two columns: the first is the prefix of fastq files (without _1.fastq.gz and _2.fastq.gz) and the second an ID for the sample:
>
```bash
cat examples/list1.txt
sample3_S1_L001 sample3_S1_L001
sample4_S1_L001 sample4_S1_L001
```
Notice that "sample3_S1_L001" appeared once. So, do not duplicate the names because the sample is paired-end (sample3_S1_L001_1.fastq.gz and sample3_S1_L001_2.fastq.gz)

Considering that sample3 and sample4 are two runs of the same animal, the sample list would be like that:

```bash
cat examples/list2.txt
sample3_S1_L001 ANIMAL1
sample4_S1_L001 ANIMAL1
sample5_S1_L001 ANIMAL2
sample6_S1_L001 ANIMAL3
```
with ANIMAL1 from sample3 and sample4. You can choose any name as IDs for the second column.

> **NOTE3**: If you do not specify the sample list (LIST_SAMPLES=nome), IVDP will create a list with all the samples in the folder and the variant calling programs will treat each sample as an independent individual.

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

#### Example 4
It mimics when you receive batches of fastq files at different times. Using the gvcf approach (gatk4GenomicsDBImport or gatk4CombineGVCFs) , the first 2 samples (examples/pe_wgs1 folder) and stop the analysis until STEP_VAR_CALL=yes to create the gvcfs:
```bash
# Run IVDP on a local server
./ivdp.sh -p params/pe_ex4a.txt
# Run IVDP on HPCC with slurm
./ivdp.sh -p params/pe_ex4a.txt -c params/configSlurm.txt
```
After that, run the other 4 samples until STEP_VAR_CALL=yes:
```bash
# Run IVDP on a local server
./ivdp.sh -p params/pe_ex4b.txt
# Run IVDP on HPCC with slurm
./ivdp.sh -p params/pe_ex4b.txt -c params/configSlurm.txt
```
When you got all the gvcfs, deactivate all the steps before STEP_VAR_CALL and activate STEP_GVCF_TO_VCF=yes to join the gvcfs and optionally you can activate STEP_VCF_FILTER=yes to filter out the final raw vcf file.
```bash
# Run IVDP on a local server
./ivdp.sh -p params/pe_ex4c.txt
# Run IVDP on HPCC with slurm
./ivdp.sh -p params/pe_ex4c.txt -c params/configSlurm.txt
```
  >**NOTE1**: Notice that  OUTPUT_NAME is the same in pe_ex4a.txt, pe_ex4b.txt, and pe_ex4c.txt. In this example, if you change the name of the analysis in any parameter files, IVDP will output the gvcfs in different folders and it will be impossible to combine the gvcfs automatically with IVDP.
 
> **NOTE2**: gatk4 genomicsDBImport [cannot accept multiple GVCFs for the same](https://gatk.broadinstitute.org/hc/en-us/articles/360035889971--How-to-Consolidate-GVCFs-for-joint-calling-with-GenotypeGVCFs/#). If you have multiple GVCFs for the same sample use gatk4CombineGVCFs like in [example 3b](#example-3)).

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

#### Example 5
WGS using two aligners (**bwa2** and **bowtie2**), three variant callers (**bcftools**, **freebayes-parallel**, and **gatk4GenomicsDBImport**). 

After filtering, the raw vcfs will be combined using 2 different strategies(**COMBINE_VCF=partial,full**). 

```bash
# Run IVDP on a local server
./ivdp.sh -p params/pe_ex5.txt
# Run IVDP on HPCC with slurm
./ivdp.sh -p params/pe_ex5.txt -c params/configSlurm.txt
```

The raw vcfs from bwa2_bcftools, bwa2_freebayes-parallel, and bwa2_gatk4GenomicsDBImport will be combined in a single vcf and the raw vcfs from bowtie2_bcftools, bowtie2_freebayes-parallel, and bowtie2_gatk4GenomicsDBImport will be combined in another vcf. The logic of COMBINE_VCF is described bellow:


![alt text](https://github.com/rodrigopsav/IVDP/blob/main/program/01.general_steps/image/example5.png)

**COMBINE_VCF=partial** keeps 
positions that appeared at least in 2 different vcfs nested to a aligner

**COMBINE_VCF=full** keeps 
positions that appeared in all vcfs nested to a aligner

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

### Example exclusive for HPCC with slurm (example 6)

#### Example 6
Same as example 2, but splitting each chromosome in chuncks of 200,000 base pairs (**BP_BY_CHROM=200000**). It will speed up the variant calling step, breaking each chromosome in chunks.

```bash
# Run IVDP on HPCC with slurm
./ivdp.sh -p params/pe_ex6.txt -c params/configSlurm.txt
```

> **NOTE1**: For real data, BP_BY_CHROM between 20000000 and 50000000 works well.
> **NOTE2**: BP_BY_CHROM does not work in local servers.

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>


### Paired-end RNAseq examples (examples 7 to 9)
The next 3 examples will use paired-end RNAseq data (IVDP/examples/pe_rnaseq). All the parameters all the same as the previous examples, with an extra step variable: **STEP_GENECOUNT=yes**. We also removed the variable **STEP_GVCF_TO_VCF** because the gvcf strategie is not used for RNAseq data.

#### Example 7
RNAseq data analysis (**ANALYSIS=2**), with parameters to include the gene counts: **STEP_GENECOUNT=yes**, **GENECOUNT_PROGRAM=featurecounts**, and **FEATURE_TYPE=exon**. We are going to use **ALIGNER_PROGRAM=hisat2** because hisat2, STAR, or gsnap are more suitable for RNAseq data (because they can handle with splice sites) compared with bowtie, bowtie2, bwa, and bwa2

```bash
# Run IVDP on a local server
./ivdp.sh -p params/pe_ex7.txt
# Run IVDP on HPCC with slurm
./ivdp.sh -p params/pe_ex7.txt -c params/configSlurm.txt
```

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

#### Example 8
Same as example 7, but using STAR aligner instead on hisat2 (**ALIGNER_PROGRAM=star**) and htseq instead of featurecounts for gene counts (**GENECOUNT_PROGRAM=featurecounts**)

```bash
# Run IVDP on a local server
./ivdp.sh -p params/pe_ex8.txt
# Run IVDP on HPCC with slurm
./ivdp.sh -p params/pe_ex8.txt -c params/configSlurm.txt
```

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

#### Example 9
Same as example 7, but using htseq and featurecounts at the same time, to get the gene count tables from both programs (**GENECOUNT_PROGRAM=htseq, featurecounts**).

```bash
# Run IVDP on a local server
./ivdp.sh -p params/pe_ex9.txt
# Run IVDP on HPCC with slurm
./ivdp.sh -p params/pe_ex9.txt -c params/configSlurm.txt
```

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

### Single-end RNAseq examples (examples 10 to 13)
The next 3 examples will use single-end RNAseq data (IVDP/examples/se_rnaseq). All the parameters all the same as examples 7 to 9, but instead using:

**EXTENSION_PE1=_1.fastq.gz**
**EXTENSION_PE2=_2.fastq.gz**

we are going to use:

**EXTENSION_SE=.fastq.gz**

Go to IVDP/examples/se_rnaseq to check that the prefix of single-end fastq files is ".fastq.gz".

We also changed the read length variables to allow smaller reads with is compatible with RNAseq technologies:

**MIN_READ_LENGTH=36**
**MAX_READ_LENGTH=70**

#### Example 10
Same as Example 7.

```bash
# Run IVDP on a local server
./ivdp.sh -p params/se_ex10.txt
# Run IVDP on HPCC with slurm
./ivdp.sh -p params/se_ex10.txt -c params/configSlurm.txt
```

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

#### Example 11
Same as Example 7 and 10, but including platypus as variant caller.

```bash
# Run IVDP on a local server
./ivdp.sh -p params/se_ex11.txt
# Run IVDP on HPCC with slurm
./ivdp.sh -p params/se_ex11.txt -c params/configSlurm.txt
```

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

#### Example 12
Same as Example 11, but we are going to get **ONLY** the gene counts table and **NOT** the vcf file, using these parameters:

STEP_MDUP=no
STEP_BQSR=no
STEP_VAR_CALL=no
STEP_VCF_FILTER=no

```bash
# Run IVDP on a local server
./ivdp.sh -p params/se_ex12.txt
# Run IVDP on HPCC with slurm
./ivdp.sh -p params/se_ex12.txt -c params/configSlurm.txt
```

> **NOTE1**:  Once STEP_VCF_FILTER=no, it does not matter the options used with **CALLER_PROGRAM** variable. IVDP will  ignore it.

> **NOTE2**: Once STEP_VCF_FILTER=no, it does not matter the options used with **MIN_DEPTH**, **MAF**, **MISSING**, and **COMBINE_VCF**. IVDP will  ignore them.

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

#### Example 13
Sample as previous example, but with 2 aligners (**ALIGNER_PROGRAM=hisat2, star**), 3 variant callers (**CALLER_PROGRAM=bcftools, freebayes-parallel, gatk4**), 2 programs for gene count (**GENECOUNT_PROGRAM=htseq, featurecounts**), with 2 strategies to combine the filtered vcfs (**COMBINE_VCF=partial, full**).

```bash
# Run IVDP on a local server
./ivdp.sh -p params/se_ex13.txt
# Run IVDP on HPCC with slurm
./ivdp.sh -p params/se_ex13.txt -c params/configSlurm.txt
```

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

## Caveats
* **Running IVDP**   
-Call IVDP from its folder. Avoid run it out of IVDP folder.   
-Set up carefully all directories path in ivdp parameter files (INPUT_DIR, OUTPUT_DIR, LIST_SAMPLES, REFERENCE, VARIANTS, ANNOTATION variables). Use absolute path to indicate them or relative path if the data and the output are inside IVDP folder (remember enter in IVDP folder before run it).   

* **Reference genome**   
-If you are running any aligner for the first time with a new reference genome, wait until the index reference genome step finishes before submit another IVDP analysis.   
E.g. you are going to start a analysis with bowtie2 using the mice genome and you had never used the combination of bowtie2 with this genome. 
Wait until the index step finishes. After that, you can submit more samples with another analysis.

* **IVDP parameter file**   
-Keep all the vertical bar "|" before variables in IVDP parameter file;   
-BP_BY_CHROM works only on HPCC with slurm;   
-CALL_BY_CHROM, THREADS and BATCH works only with IVDP in a local server.   
-Do not use **DOWNSAMPLING_FASTQ** and **DOWNSAMPLING_BAM** at the same time. Choose only one.  
-Do not use freebayes-parallel on HPCC with slurm if BP_BY_CHROM is different of "all".   
-Use gatk4CombineGVCFs only if:   
  * gatk4GenomicsDBImport fails;   
  * you have repeated samples from the same individual (i.e. various fastq files for same individual). Go to [Example 3](https://stackedit.io/app#example-3) to understand how to handle with repeated runs.   

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

## Troubleshooting

### 1. Installation
-If trimmomatic, picard or gatk failled, try to install manually open java:
```
conda install -y -n ivdp -c anaconda openjdk 
```

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

### 2. Problems to run IVDP
```diff
- ERROR: conda ivdp environment was not found.
```
First check if ivdp, ivdp2 or r-env4.0 are activated. If so, use:
```
conda deactivate
```
Otherwise, [Install IVDP dependencies](#install-ivdp-dependencies).

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

### 3. My analysis failed in a certain step
Let's pretend that a analysis failed on mark duplicated reads step. Set up all steps before mdup equal to "no", like that:

|STEP_QC_TRIM=**no**   
|STEP_ALIGNMENT=**no**   
|STEP_MDUP=yes   
|STEP_BQSR=yes   
|STEP_VAR_CALL=yes   
|STEP_GVCF_TO_VCF=yes   
|STEP_VCF_FILTER=yes   

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>

### 4. Some samples failed in a certain step
Let's suppose you are running [Example 5](#example-5) with 10 samples and 2 (samples7 and sample9) of them failed at mark duplicated reads step. Instead of running the analysis for all the samples again, you can specify a list of samples that you want to run:

```bash
vim IVDP/examples/list3.txt
sample7 sample7
sample9 sample9
```

**-Run only this samples for the failed step:**

|LIST_SAMPLES=**~/IVDP/examples/list3.txt**   

|STEP_QC_TRIM=**no**   
|STEP_ALIGNMENT=**no**   
|STEP_MDUP=**yes**   
|STEP_BQSR=**no**   
|STEP_VAR_CALL=**no**   
|STEP_GVCF_TO_VCF=**no**   
|STEP_VCF_FILTER=**no**   

**-After finish, run the remaining steps, without any list of samples:**

|LIST_SAMPLES=**none**   

|STEP_QC_TRIM=no   
|STEP_ALIGNMENT=no   
|STEP_MDUP=no   
|STEP_BQSR=**yes**   
|STEP_VAR_CALL=**yes**   
|STEP_GVCF_TO_VCF=**yes**   
|STEP_VCF_FILTER=**yes**   

<div align="right">
    <b><a href="#ivdp-integrated-variant-discovery-pipeline">↥ back to top</a></b>
</div>
