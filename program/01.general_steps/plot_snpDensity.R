##############################################
##### PLOT SNP DENSITY ALONG CHROMOSOMES #####
##############################################
my_packages <- c("data.table", "ggplot2", "gridExtra", "rMVP") # Specify your packages
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])] # Extract not installed packages
if(length(not_installed)) install.packages(not_installed, repos = "https://cran.microsoft.com/", dependencies=TRUE)

##### Set initial parameters #####
library(data.table)
library(ggplot2)
library(gridExtra)

inputFile=inputFile
outputPath=outputPath
aligner=aligner
caller=caller
window <- 500000 # window length in bp
setwd(outputPath)

#---------- VARIANT POSITIONS ----------#
variant_positions <- fread(inputFile, header=F, sep = "\t")
variant_positions <- as.data.frame(variant_positions)
colnames(variant_positions) <- c("Chromosome", "Position", "ID")

# Renaming positions without rsID
variant_positions[,3] <-ifelse(variant_positions[,3]==".", 
                               paste0(variant_positions[,1], ":", variant_positions[,2]), 
                               variant_positions[,3])

# Change column order
variant_positions <- as.data.frame(variant_positions[ , c(3,1,2)])

# Convert character to factor
variant_positions[ ,1] <- as.factor(variant_positions[ ,1])
variant_positions[ ,2] <- as.factor(variant_positions[ ,2])

#---------- Plot SNP density ----------#
### https://www.rdocumentation.org/packages/rMVP/versions/0.99.14.1
#install.packages("bigmemory.sri")
#install.packages("rMVP")
suppressMessages(library(rMVP))

# All the genome
print("Creating plots: snp density")
MVP.Report(variant_positions[, c(1:3)], plot.type="d", 
           col=c("darkgreen", "yellow", "red"),
           bin.size = window,
           file.output = TRUE, outpath = outputPath,
           memo = paste0("plot_DensitySNPs_", aligner, "_", caller), 
           file.type = "pdf", dpi=300, width=12)

# Chromosome 23
# MVP.Report(variant_positions[variant_positions(vcf1[ ,2] == "23"), c(1:3)], plot.type="d", 
#            col=c("darkgreen", "yellow", "red"),
#            bin.size = window,
#            file.output = TRUE, outpath = outputPath,
#            memo = paste0("plot_DensitySNPs_chrom23_", aligner, "_", caller), 
#            file.type = "pdf", dpi=300)
