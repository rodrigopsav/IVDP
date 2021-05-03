#!/bin/bash

################################
##### INSTALL DEPENDENCIES #####
################################

##### SET PARAMETERS
usage() { echo "Usage: $0 -d path/to/install/dependencies/programs
       -d: directory path to install dependencies programs
       " 1>&2; exit 1; }

while getopts :d: option; do
   case "${option}" in
   d) INSTALL_FOLDER=${OPTARG};;
   *) usage;;
   esac
done
#shift "$((OPTIND-1))"

##### TEST INSTALL_FOLDER VARIABLE
if [[ -z "$INSTALL_FOLDER" ]]; then
   echo "Error: -d flag is empty"
   usage
   exit 1

else

   export INSTALL_FOLDER=$(readlink -f $INSTALL_FOLDER)
   
   if [[ ! -d "$INSTALL_FOLDER" ]]; then
      echo "ERROR: wrong directory path. Please check -d flag"
      echo "Aborting analysis"
      usage
      exit 1
   fi
fi
wait

echo "Install programs in: "$INSTALL_FOLDER


##### Make Installation folder
#export INSTALL_FOLDER=$INSTALL_FOLDER
#mkdir -p $INSTALL_FOLDER


echo "
#--------------------------------#
##### INSTALLING MINICONDA 3 #####
#--------------------------------#
"

if ! command -v conda &> /dev/null; then
   echo "conda CANNOT BE FOUND"
   
   # Make Installation folder
   export INSTALL_FOLDER=$INSTALL_FOLDER
   mkdir -p $INSTALL_FOLDER

   cd $INSTALL_FOLDER
   if [[ "$OSTYPE" == "linux-gnu"* ]]; then
      wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
      chmod +x Miniconda3-latest-Linux-x86_64.sh
      # https://docs.anaconda.com/anaconda/install/silent-mode/
      ./Miniconda3-latest-Linux-x86_64.sh -b -p $INSTALL_FOLDER/miniconda3 -f
      source $INSTALL_FOLDER/miniconda3/bin/activate
      conda init bash
      #https://www.kangzhiq.com/2020/05/02/how-to-activate-a-conda-environment-in-a-remote-machine/
      eval "$(conda shell.bash hook)"
      rm Miniconda3-latest-Linux-x86_64.sh
   elif [[ "$OSTYPE" == "darwin"* ]]; then
      wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
      chmod +x Miniconda3-latest-MacOSX-x86_64.sh
      ./Miniconda3-latest-MacOSX-x86_64.sh -b -p $INSTALL_FOLDER/miniconda3 -f
      source $INSTALL_FOLDER/miniconda3/bin/activate
      conda init bash
      #https://www.kangzhiq.com/2020/05/02/how-to-activate-a-conda-environment-in-a-remote-machine/
      eval "$(conda shell.bash hook)"
      rm Miniconda3-latest-MacOSX-x86_64.sh
   else
      "Error: This is neither linux nor MacOS system"
       exit 1
   fi

else
   echo "conda is installed" 
fi
   
### Configure conda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels anaconda
conda config --add channels conda-forge

### Create conda environment to save Bioinformatic programs
#https://www.kangzhiq.com/2020/05/02/how-to-activate-a-conda-environment-in-a-remote-machine/
conda init bash
eval "$(conda shell.bash hook)"
conda create -y --name ivdp python=3.8
conda create -y --name ivdp2 python=2.7
conda create -y --name r-env4.0 -c conda-forge r-base=4.0.3 python=3.8

conda activate ivdp

# NOTE: If the conda ivdp environment is not created, some programs, like fastqc, won't work 
# due to conflicts with the base conda environment 


echo "
#-----------------------------#
##### INSTALLING PARALLEL #####
#-----------------------------#
"
conda install -y -n ivdp -c conda-forge parallel=20200722


echo "
#-----------------------------#
##### INSTALLING DOS2UNIX #####
#-----------------------------#
"
conda install -y -n ivdp -c trent dos2unix=7.4.1


echo "
#---------------------------#
##### INSTALLING RENAME #####
#---------------------------#
"
conda install -y -n ivdp -c bioconda rename=1.601


echo "
#------------------------#
##### INSTALLING SED #####
#------------------------#
"
conda install -y -n ivdp -c conda-forge sed=4.8


echo "
#------------------------#
##### INSTALLING GIT #####
#------------------------#
"
conda install -y -n ivdp -c anaconda git=2.28.0


echo "
#-------------------------#
##### INSTALLING PIGZ #####
#-------------------------#
"
conda install -y -n ivdp -c anaconda pigz=2.4


echo "
#-----------------------------#
##### INSTALLING DATAMASH #####
#-----------------------------#
"
conda install -y -n ivdp -c anaconda datamash=1.7


echo "
#-------------------------------#
##### INSTALLING SRATOOLKIT #####
#-------------------------------#
"
conda install -y -n ivdp -c bioconda sra-tools=2.10.8


echo "
#----------------------------------#
##### INSTALLING ENTREZ-DIRECT #####
#----------------------------------#
"
conda install -y -n ivdp -c bioconda entrez-direct=13.9


echo "
#---------------------------#
##### INSTALLING FASTQC #####
#---------------------------#
"
conda install -y -n ivdp -c bioconda fastqc=0.11.9


echo "
#--------------------------------#
##### INSTALLING FASTQ_UTILS #####
#--------------------------------#
"
conda install -y -n ivdp -c bioconda fastq_utils=0.24.1


echo "
#--------------------------#
##### INSTALLING SEQTK #####
#--------------------------#
"
conda install -y -n ivdp -c bioconda seqtk=1.3


echo "
#----------------------------#
##### INSTALLING MULTIQC #####
#----------------------------#
"
conda install -y -n ivdp -c bioconda -c conda-forge multiqc=1.9


echo "
#---------------------------------#
##### INSTALLING TRIMMOMATICS #####
#---------------------------------#
"
conda install -y -n ivdp -c bioconda trimmomatic=0.39


echo "
#--------------------------#
##### INSTALLING BBMAP #####
#--------------------------#
"
conda install -y -n ivdp -c bioconda bbmap=38.90


echo "
#----------------------------#
##### INSTALLING BOWTIE ######
#----------------------------#
"
conda install -y -n ivdp -c bioconda bowtie=1.3.0


echo "
#----------------------------#
##### INSTALLING BOWTIE2 #####
#----------------------------#
"
conda install -y -n ivdp -c bioconda bowtie2=2.4.2


echo "
#------------------------#
##### INSTALLING BWA #####
#------------------------#
"
conda install -y -n ivdp -c bioconda bwa=0.7.17


echo "
#-----------------------------#
##### INSTALLING BWA-MEM2 #####
#-----------------------------#
"
conda install -y -n ivdp -c bioconda bwa-mem2=2.2.1


echo "
#-------------------------------#
##### INSTALLING GMAP-GSNAP #####
#-------------------------------#
"
conda install -y -n ivdp -c anaconda zlib=1.2.11
conda install -y -n ivdp -c bioconda gmap=2021.03.08


echo "
#---------------------------#
##### INSTALLING HISAT2 #####
#---------------------------#
"
#conda install -y -n ivdp -c biobuilds hisat2=2.2.1
conda install -y -n ivdp -c bioconda hisat2=2.2.1


echo "
#-------------------------#
##### INSTALLING STAR #####
#-------------------------#
"
conda install -y -n ivdp -c bioconda star=2.7.8a


echo "
#-----------------------------------------------#
##### INSTALLING SAMTOOLS, BCFTOOLS, HTSLIB #####
#-----------------------------------------------#
"
conda install -y -n ivdp -c bioconda htslib=1.9
conda install -y -n ivdp2 -c bioconda htslib=1.9
conda install -y -n ivdp -c bioconda samtools=1.9
conda install -y -n ivdp -c bioconda bcftools=1.9


echo "
#-----------------------------#
##### INSTALLING SAMBAMBA #####
#-----------------------------#
"
conda install -y -n ivdp -c bioconda sambamba=0.7.1


echo "
#---------------------------#
##### INSTALLING PICARD #####
#---------------------------#
"
conda install -y -n ivdp -c bioconda picard=2.25.2


echo "
#-----------------------------------#
##### INSTALLING GATK VERSION 4 #####
#-----------------------------------#
"
conda install -y -n ivdp -c bioconda gatk4=4.2.0.0


echo "
#------------------------------#
##### INSTALLING FREEBAYES #####
#------------------------------#
"
conda install -y -n ivdp2 -c bioconda freebayes=1.3.2
#./conda config --add channels defaults
#./conda config --add channels bioconda
#./conda config --add channels conda-forge
#./conda install freebayes=1.3.2-0


echo "
#-----------------------------#
##### INSTALLING PLATYPUS #####
#-----------------------------#
"
conda install -y -n ivdp2 -c bioconda platypus-variant=0.8.1.2


echo "
#-----------------------------#
##### INSTALLING VARSCAN2 #####
#-----------------------------#
"
conda install -y -n ivdp -c bioconda varscan=2.4.4
#cd $INSTALL_FOLDER
#git clone https://github.com/dkoboldt/varscan.git


echo "
#-------------------------#
##### INSTALLING PERL #####
#-------------------------#
"
#conda install -y -n ivdp -c anaconda perl


echo "
#-----------------------------#
##### INSTALLING VCFTOOLS #####
#-----------------------------#
"
conda install -y -n ivdp -c bioconda vcftools=0.1.16
conda install -y -n ivdp -c bioconda perl-vcftools-vcf


echo "
#--------------------------#
##### INSTALLING plink #####
#--------------------------#
"
conda install -y -n ivdp -c bioconda plink=1.90b6.18
conda install -y -n ivdp -c bioconda plink2=v2.00a2.3LM


echo "
#-----------------------------#
##### INSTALLING BEDTOOLS #####
#-----------------------------#
"
conda install -y -n ivdp -c bioconda bedtools=2.29.2


echo "
#--------------------------#
##### INSTALLING HTSEQ #####
#--------------------------#
"
#conda install -y -n ivdp -c anaconda numpy
#conda install -y -n ivdp -c bioconda pysam 
#conda install -y -n ivdp -c bioconda htseq
mkdir -p $INSTALL_FOLDER/tmp
TMPDIR=$INSTALL_FOLDER/tmp pip install htseq==0.12.4
wait
rm -r $INSTALL_FOLDER/tmp


echo "
#----------------------------------#
##### INSTALLING FEATURECOUNTS #####
#----------------------------------#
"
conda install -y -n ivdp -c bioconda subread=2.0.1


echo "
#---------------------------------------#
##### INSTALLING SNPEFF AND SNPSIFT #####
#---------------------------------------#
"
conda install -y -n ivdp -c bioconda snpeff=5.0 
conda install -y -n ivdp -c bioconda snpsift=4.3.1t


echo "
#---------------------------#
##### INSTALLING VCFLIB #####
#---------------------------#
"
conda install -y -n ivdp -c bioconda vcflib=1.0.0_rc2
conda install -y -n ivdp2 -c bioconda vcflib=1.0.0_rc2


echo "
#--------------------------------#
##### INSTALLING ghostscript #####
#--------------------------------#
"
conda install -y -n ivdp -c conda-forge ghostscript


echo "
#------------------------------#
##### INSTALLING R Program #####
#------------------------------#
"
conda install -y -n r-env4.0 "libblas=*=*openblas"
conda install -y -n r-env4.0 -c conda-forge r-dplyr
conda install -y -n r-env4.0 -c conda-forge r-data.table
conda install -y -n r-env4.0 -c conda-forge r-ggplot2
conda install -y -n r-env4.0 -c conda-forge r-gridextra
conda install -y -n r-env4.0 -c conda-forge r-bigmemory.sri
conda install -y -n r-env4.0 -c conda-forge r-rmvp
conda install -y -n r-env4.0 -c conda-forge r-formattable
conda install -y -n r-env4.0 -c conda-forge r-htmltools
conda install -y -n r-env4.0 -c conda-forge r-webshot


echo "
#----------------------------#
##### INSTALLING FIREFOX #####
#----------------------------#
"
conda install -y -n ivdp -c conda-forge firefox=87.0


echo "
#----------------------------#
##### INSTALLING OPENJDK #####
#----------------------------#
#"
conda install -y -n ivdp -c anaconda openjdk=8.0.282


##### Deactivate ivdp conda environment
conda deactivate


echo "
#----------------------------------#
##### CHECK INSTALLED PROGRAMS #####
#----------------------------------#
"

echo "#############################################################################"
echo "############## CHECK IF ALL THE PROGRAMS ARE INSTALLED ######################"
echo "If any installed program failed, re-install a previous version. To do that:"
echo "Example using fastqc package"
echo "Search for versions: conda search fastqc --info"
echo "Choose a previous version: eg. 0.11.8"
echo "Install: conda install --force-reinstall -y -n ivdp -c bioconda fastqc=0.11.8"
echo "#############################################################################"
echo


### Activate conda ivdp environment 
eval "$(conda shell.bash hook)"
conda activate ivdp
conda activate --stack ivdp2
conda activate --stack r-env4.0


###
if ! command -v parallel &> /dev/null; then
   echo "parallel CANNOT BE FOUND"
else
   echo "parallel is installed" 
fi

###
if ! command -v dos2unix &> /dev/null; then
   echo "dos2unix CANNOT BE FOUND"
else
   echo "dos2unix is installed" 
fi

###
if ! command -v rename &> /dev/null; then
   echo "rename CANNOT BE FOUND"
else
   echo "rename is installed" 
fi

###
if ! command -v sed &> /dev/null; then
   echo "sed CANNOT BE FOUND"
else
   echo "sed is installed" 
fi

###
if ! command -v git &> /dev/null; then
   echo "git CANNOT BE FOUND"
else
   echo "git is installed" 
fi

###
if ! command -v pigz &> /dev/null; then
   echo "pigz CANNOT BE FOUND"
else
   echo "pigz is installed" 
fi

###
if ! command -v datamash &> /dev/null; then
   echo "datamash CANNOT BE FOUND"
else
   echo "datamash is installed" 
fi

###
if ! command -v prefetch &> /dev/null; then
   echo "sra-toolkit CANNOT BE FOUND"
else
   echo "sra-toolkit is installed" 
fi

###
if ! command -v esearch &> /dev/null; then
   echo "entrez-direct CANNOT BE FOUND"
else
   echo "entrez-direct is installed" 
fi

###
if ! command -v fastqc &> /dev/null; then
   echo "fastqc CANNOT BE FOUND"
else
   echo "fastqc is installed" 
fi

###
if ! command -v fastq_info &> /dev/null; then
   echo "fastq_utils CANNOT BE FOUND"
else
   echo "fastq_utils is installed" 
fi

###
if ! command -v seqtk &> /dev/null; then
   echo "seqtk CANNOT BE FOUND"
else
   echo "seqtk is installed" 
fi

###
if ! command -v multiqc &> /dev/null; then
   echo "multiqc CANNOT BE FOUND"
else
   echo "multiqc is installed" 
fi

###
if ! command -v trimmomatic &> /dev/null; then
   echo "trimmomatic CANNOT BE FOUND"
else
   echo "trimmomatic is installed" 
fi

###
if ! command -v bbmap.sh &> /dev/null; then
   echo "bbmap.sh CANNOT BE FOUND"
else
   echo "bbmap.sh is installed" 
fi

###
if ! command -v bowtie &> /dev/null; then
   echo "bowtie CANNOT BE FOUND"
else
   echo "bowtie is installed" 
fi

###
if ! command -v bowtie2 &> /dev/null; then
   echo "bowtie2 CANNOT BE FOUND"
else
   echo "bowtie2 is installed" 
fi

###
if ! command -v bwa &> /dev/null; then
   echo "bwa CANNOT BE FOUND"
else
   echo "bwa is installed" 
fi

###
if ! command -v bwa-mem2 &> /dev/null; then
   echo "bwa-mem2 CANNOT BE FOUND"
else
   echo "bwa-mem2 is installed" 
fi

###
if ! command -v gsnap &> /dev/null; then
   echo "gsnap CANNOT BE FOUND"
else
   echo "gsnap is installed" 
fi

###
if ! command -v hisat2 &> /dev/null; then
   echo "hisat2 CANNOT BE FOUND"
else
   echo "hisat2 is installed" 
fi

###
if ! command -v STAR &> /dev/null; then
   echo "STAR CANNOT BE FOUND"
else
   echo "STAR is installed" 
fi

###
if ! command -v sambamba &> /dev/null; then
   echo "sambamba CANNOT BE FOUND"
else
   echo "sambamba is installed" 
fi

###
if ! command -v samtools &> /dev/null; then
   echo "samtools CANNOT BE FOUND"
else
   echo "samtools is installed" 
fi

###
if ! command -v bgzip &> /dev/null; then
   echo "htslib CANNOT BE FOUND"
else
   echo "htslib is installed" 
fi

###
if ! command -v picard &> /dev/null; then
   echo "picard CANNOT BE FOUND"
else
   echo "picard is installed" 
fi

###
if ! command -v bcftools &> /dev/null; then
   echo "bcftools CANNOT BE FOUND"
else
   echo "bcftools is installed" 
fi

###
if ! command -v gatk &> /dev/null; then
   echo "gatk CANNOT BE FOUND"
else
   echo "gatk is installed" 
fi

###
if ! command -v varscan &> /dev/null; then
   echo "varscan CANNOT BE FOUND"
else
   echo "varscan is installed" 
fi

###
if ! command -v freebayes &> /dev/null; then
   echo "freebayes CANNOT BE FOUND"
else
   echo "freebayes is installed" 
fi

###
if ! command -v freebayes-parallel &> /dev/null; then
   echo "freebayes-parallel CANNOT BE FOUND"
else
   echo "freebayes-parallel is installed" 
fi

###
if ! command -v platypus &> /dev/null; then
   echo "platypus CANNOT BE FOUND"
else
   echo "platypus is installed" 
fi

###
if ! command -v vcftools &> /dev/null; then
   echo "vcftools CANNOT BE FOUND"
else
   echo "vcftools is installed" 
fi

###
if ! command -v plink &> /dev/null; then
   echo "plink COULD NOT BE FOUND"
else
   echo "plink is installed" 
fi

###
if ! command -v plink2 &> /dev/null; then
   echo "plink2 COULD NOT BE FOUND"
else
   echo "plink2 is installed" 
fi

###
if ! command -v bedtools &> /dev/null; then
   echo "bedtools CANNOT BE FOUND"
else
   echo "bedtools is installed" 
fi

###
if ! command -v htseq-count &> /dev/null; then
   echo "htseq-count CANNOT BE FOUND"
else
   echo "htseq-count is installed" 
fi

###
if ! command -v featureCounts &> /dev/null; then
   echo "featureCounts CANNOT BE FOUND"
else
   echo "featureCounts is installed" 
fi

###
if ! command -v snpEff &> /dev/null; then
   echo "snpEff CANNOT BE FOUND"
else
   echo "snpEff is installed" 
fi

###
if ! command -v SnpSift &> /dev/null; then
   echo "SnpSift CANNOT BE FOUND"
else
   echo "SnpSift is installed" 
fi

###
if ! command -v vcffilter &> /dev/null; then
   echo "vcflib CANNOT BE FOUND"
else
   echo "vcflib is installed" 
fi

###
if ! command -v gs &> /dev/null; then
   echo "ghostscript COULD NOT BE FOUND"
else
   echo "ghostscript is installed" 
fi

###
if ! command -v R &> /dev/null; then
   echo "R Program COULD NOT BE FOUND"
else
   echo "R Program is installed" 
fi


##### Deactivate ivdp conda environment
conda deactivate
conda deactivate
conda deactivate
echo


echo "#################################################################################"
echo "# TO RUN THESE PROGRAMS FOR OTHER APPLICATION, ACTIVATE IVDP CONDA ENVIRONMENTS #"
echo "############ TO ACTIVATE, TYPE ON TERMINAL: conda activate ivdp #################"
echo "############ TO ACTIVATE, TYPE ON TERMINAL: conda activate ivdp2 ################"
echo "#################################################################################"

echo "
#--------------------------------#
##### INSTALLATION COMPLETED #####
#--------------------------------#
"

