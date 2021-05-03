#########################################
########## VCF STATISTICS PLOT ##########
#########################################
my_packages <- c("data.table", "ggplot2", "gridExtra") # Specify your packages
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])] # Extract not installed packages
if(length(not_installed)) install.packages(not_installed)

library(data.table)
library(ggplot2)
library(gridExtra)

inputPath=inputPath
aligner=aligner
caller=caller
step=step
chromNames=chromNames
varCountData=file1
freqData=file2
vmissData=file3
smissData=file4
hetData=file5
setwd(inputPath)

#------------------ VARIANT COUNTS ----------------#
chromNames <- fread(chromNames, header=F)
chromNames <- unlist(chromNames)
counts <- fread(varCountData, header = F)
colnames(counts) <- c("chrom", "variant", "n")

print("Creating plots: variant counts")
p1 <- ggplot(counts, aes(x=chrom, y=n, fill=variant)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Variants by chromosome", subtitle = bquote(bold("Aligner: ")~.(aligner)~~~~bold("Caller: ")~.(caller)~~~~bold("File: ")~.(step))) + 
  scale_x_continuous(labels = chromNames, breaks = seq(1, length(chromNames), by=1)) + 
  scale_fill_manual(breaks =c("snp", "indel"), values=c("lightblue3", "pink3"), name="Variant\ntype") +
  theme_classic() + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        panel.border = element_blank(), 
        axis.text.x = element_text(size=14, angle = 0), 
        axis.text.y = element_text(size=14), 
        axis.title.x = element_text(size=14, face="bold"), 
        axis.title.y = element_text(size=14, face="bold"), 
        plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        plot.subtitle = element_text(size=15, hjust = 0.5), 
        legend.title = element_text(colour="black", size=15, face="bold"), 
        legend.text = element_text(size=14)) + 
  xlab("Chromosomes") + 
  ylab("Counts")

#---------- MINOR ALLELE FREQUENCY (MAF) ----------#
freq <- fread(freqData, header=T, select=c(5), sep = "\t")
freq <- as.data.frame(freq)
freq_alt <- as.numeric(unlist(freq))
maf <- ifelse(freq_alt <= 0.5, freq_alt, 1 - freq_alt)

# hist(maf, col = "skyblue3",
#      main="MAF",
#      prob = F,
#      xlab = "MAF",
#      breaks = c(seq(0, 0.5, 0.01)))


print("Creating plots: allele frequencies")
p2 <- ggplot(data.frame(maf=maf), aes(x=maf)) + 
  geom_histogram(aes(y=..density..), breaks = c(seq(0, 0.5, 0.01)), na.rm = T, fill="skyblue3", color="skyblue3", alpha=0.9) +
  ggtitle("Minor allele frequency", subtitle = bquote(bold("Aligner: ")~.(aligner)~~~~bold("Caller: ")~.(caller)~~~~bold("File: ")~.(step))) +
  theme_light() +
  theme(panel.grid.minor.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_line(colour = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        plot.title = element_text(size=12, hjust = 0.5, face="bold"),
        plot.subtitle = element_text(size=10, hjust = 0.5),
        legend.position = "none") +
  xlab("MAF") + 
  ylab("Frequency")

# https://www.datanovia.com/en/lessons/ggplot-ecdf/
p3 <- ggplot(data.frame(maf=maf), aes(x=maf)) +
  stat_ecdf(aes(x=maf),geom = "step", pad = FALSE, na.rm = T, size = 0.5, color="skyblue3") +
  scale_x_continuous(breaks = c(seq(0, 0.5, 0.05))) +
  scale_y_continuous(breaks = c(seq(0, 1, 0.1))) +
  ggtitle("Cumulative Density Function MAF", subtitle = bquote(bold("Aligner: ")~.(aligner)~~~~bold("Caller: ")~.(caller)~~~~bold("File: ")~.(step))) +
  theme_classic() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_line(colour = "black"),
        axis.text.x = element_text(size=12, angle = 90),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        plot.title = element_text(size=12, hjust = 0.5, face="bold"),
        plot.subtitle = element_text(size=10, hjust = 0.5),
        legend.position = "none") +
  xlab("MAF") + 
  ylab("Cumulative Frequency")

#---------- GENOTYPE-BASED MISSING ----------#
vmiss <- fread(vmissData, header=T, sep = "\t")
vmiss <- as.data.frame(vmiss)
# hist(vmiss[,5], col = "salmon2",
#      main = "Genotype-based missing calls",
#      prob = F,
#      xlab = "Proportion of missing calls",
#      breaks = c(seq(0, 1, 0.05)))


print("Creating plots: missing variants by position")
p4 <- ggplot(vmiss, aes(x=F_MISS)) + 
  geom_histogram(aes(y=..density..), breaks=c(seq(0, 1, 0.05)), fill="palegreen3", color="palegreen3", alpha=0.9) +
  ggtitle("Genotype-based missing calls", subtitle = bquote(bold("Aligner: ")~.(aligner)~~~~bold("Caller: ")~.(caller)~~~~bold("File: ")~.(step))) +
  theme_light() +
  theme(panel.grid.minor.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_line(colour = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        plot.title = element_text(size=12, hjust = 0.5, face="bold"),
        plot.subtitle = element_text(size=10, hjust = 0.5),
        legend.position = "none") +
  xlab("Proportion of missing calls") + 
  ylab("Frequency")

p5 <- ggplot(vmiss, aes(x=F_MISS)) +
  stat_ecdf(aes(x=F_MISS),geom = "step", pad = FALSE, na.rm = T, size = 0.5, color="palegreen3") +
  scale_x_continuous(breaks = c(seq(0, 1, 0.1))) +
  scale_y_continuous(breaks = c(seq(0, 1, 0.1))) +
  ggtitle("CDF Genotype-based missing calls", subtitle = bquote(bold("Aligner: ")~.(aligner)~~~~bold("Caller: ")~.(caller)~~~~bold("File: ")~.(step))) +
  theme_classic() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_line(colour = "black"),
        axis.text.x = element_text(size=12, angle = 90),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        plot.title = element_text(size=12, hjust = 0.5, face="bold"),
        plot.subtitle = element_text(size=10, hjust = 0.5),
        legend.position = "none") +
  xlab("Proportion of missing calls") + 
  ylab("Cumulative Frequency")

#---------- SAMPLE-BASED MISSING ----------#
smiss <- fread(smissData, header=T, select=c(1,4), sep = "\t")
smiss <- as.data.frame(smiss)
colnames(smiss) <- c("ID", "percent")
smiss[,1] <- as.character(smiss[,1])

# hist(smiss[,2], col = "salmon2",
#      main = "Sample-based missing calls",
#      prob = F,
#      xlab = "Proportion of missing calls",
#      breaks = c(seq(0, 1, 0.05)))

print("Creating plots: missing variants by sample")
p6 <- ggplot(smiss, aes(x=percent)) + 
  geom_histogram(breaks=c(seq(0, 1, 0.05)), na.rm = T, fill="salmon2", color="salmon2", alpha=0.9) +
  ggtitle("Sample-based missing calls", subtitle = bquote(bold("Aligner: ")~.(aligner)~~~~bold("Caller: ")~.(caller)~~~~bold("File: ")~.(step))) +
  theme_light() +
  theme(panel.grid.minor.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_line(colour = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        plot.title = element_text(size=12, hjust = 0.5, face="bold"),
        plot.subtitle = element_text(size=10, hjust = 0.5),
        legend.position = "none") +
  xlab("Proportion of missing calls") + 
  ylab("Frequency")

p7 <- ggplot(smiss, aes(x = ID, y = percent)) +
  geom_segment(aes(x = ID, y = 0, xend = ID, yend = percent), color="grey") +
  geom_point( color="salmon2", size=4) +
  ggtitle("Missing calls by sample", subtitle = bquote(bold("Aligner: ")~.(aligner)~~~~bold("Caller: ")~.(caller)~~~~bold("File: ")~.(step))) +
  theme_light() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_line(colour = "black"),
        axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        plot.title = element_text(size=12, hjust = 0.5, face="bold"),
        plot.subtitle = element_text(size=10, hjust = 0.5),
        legend.position = "none") +
  xlab("") + 
  ylab("Proportion of missing calls")

#---------- INBREEDING ----------#
inbreed <- fread(hetData, header=T, select = c(2,6))
inbreed <- as.data.frame(inbreed)
colnames(inbreed) <- c("ID", "F")
inbreed[,1] <- as.character(inbreed[,1])
inbreed$index <- 1:dim(inbreed)[1]
inbreed$flag <- ifelse(inbreed$F < 0, inbreed$flag <- 'negative', inbreed$flag <- 'positive')

# hist(inbreed[,2], col = "seashell3",
#      main = "Inbreeding coefficient",
#      prob = F,
#      xlab = "Inbreeding coefficient (F)",
#      breaks = 10)

print("Creating plots: inbreeding")
p8 <- ggplot(inbreed, aes(x=F)) + 
  geom_histogram(bins=20, na.rm = T, fill="seashell3", color="seashell3", alpha=0.9) +
  ggtitle("Inbreeding coefficient (F)", subtitle = bquote(bold("Aligner: ")~.(aligner)~~~~bold("Caller: ")~.(caller)~~~~bold("File: ")~.(step))) +
  theme_light() +
  theme(panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        plot.title = element_text(size=12, hjust = 0.5, face="bold"),
        plot.subtitle = element_text(size=10, hjust = 0.5),
        legend.position = "none") +
  xlab("Inbreeding coefficient (F)") + 
  ylab("Frequency")

p9 <- ggplot(inbreed, aes(x = ID, y = F, color = flag)) +
  geom_segment(aes(x = ID, y = 0, xend = ID, yend = F), color="grey") +
  geom_point( color="seashell3", size=4) +
  geom_point() +
  ggtitle("Inbreeding coefficient (F) by individual", subtitle = bquote(bold("Aligner: ")~.(aligner)~~~~bold("Caller: ")~.(caller)~~~~bold("File: ")~.(step))) +
  theme_light() +
  theme(panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size=12, face="bold"),
    axis.title.y = element_text(size=12, face="bold"),
    plot.title = element_text(size=12, hjust = 0.5, face="bold"),
    plot.subtitle = element_text(size=10, hjust = 0.5),
    legend.position = "none") +
  xlab("") + 
  ylab("Inbreeding coefficient (F)")


#plots <- list(p1, p2, p3, p5, p4, p6)
#margin = theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
#grid.arrange(grobs=lapply(plots,"+", margin) , ncol=2)


print("Saving plots")
plots_page1 <- list(p1)
margin = theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
g1 <- arrangeGrob(grobs=lapply(plots_page1,"+", margin), ncol=1)
suppressMessages(
  ggsave(file=paste0("plot_stats_", aligner, "_", caller, "_page1.pdf"), g1,
         dpi=300,
         width = 12, height = 8)
)

plots_page2 <- list(p2, p3, p4, p5, p6, p8)
margin = theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
g2 <- arrangeGrob(grobs=lapply(plots_page2,"+", margin), ncol=2)
suppressMessages(
ggsave(file=paste0("plot_stats_", aligner, "_", caller, "_page2.pdf"), g2,
       dpi=300,
       width = 12, height = 16)
)

plots_page3 <- list(p7, p9)
margin = theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
g3 <- arrangeGrob(grobs=lapply(plots_page3,"+", margin), ncol=1)
suppressMessages(
ggsave(file=paste0("plot_stats_", aligner, "_", caller, "_page3.pdf"), g3,
       dpi=300,
       width = 40, height = 12)
)

