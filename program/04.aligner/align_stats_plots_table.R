#########################################
########## BAM STATISTICS PLOT ##########
#########################################
my_packages <- c("data.table", "bit64", "formattable", "htmltools", "webshot") # Specify your packages
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])] # Extract not installed packages
if(length(not_installed)) install.packages(not_installed, repos = "https://cran.microsoft.com/", dependencies=TRUE)

suppressMessages(library(data.table))
suppressMessages(library(bit64))
library(formattable)
library(htmltools)
library(webshot)
suppressMessages(webshot::install_phantomjs())

inputPath=inputPath
aligner=aligner
setwd(inputPath)

bam_stats <- fread(paste0(aligner, "_align_stats.txt"), header = T)
bam_stats <- as.data.frame(bam_stats)
bam_stats[,1] <- as.character(bam_stats[,1])
bam_stats[2:13] <- sapply(bam_stats[2:13],as.numeric)
bam_stats_variables <-bam_stats[ , -1]
bam_stats_variables[1:12] <- sapply(bam_stats_variables[1:12],as.numeric)

title <- c("Depth of Coverage", "Breadth of Coverage", "Depth for covered regions",
           "Total number of reads", "Number of mapped reads", "Number of unmapped reads",
           "Number of paired reads", "Number of bases mapped", "Average of read length",
           "Average quality of reads", "Percentage of mapped reads", "Percentage of unmapped reads")

ylabel <- c("Depth of Coverage (x)", "Breadth of Coverage (%)", "Depth of Coverage for covered regions (x)",
            "Total number of reads", "Number of mapped reads", "Number of unmapped reads",
            "Number of paired reads", "Number of bases mapped", "Average read length",
            "Average quality of reads", "Mapped reads (%)", "Unmapped reads (%)")

color <- c("#88CCEE", "#CC6677", "#DDCC77",
           "#117733", "#332288", "#AA4499",
           "#44AA99", "#999933", "882255",
           "#661100", "#6699CC", "#888888")

# bam_stats %>% dplyr::select(Sample, everything()) %>% 
#   tidyr::gather("id", "value",2:13) %>% 
#   ggplot(., aes(x = id, y = value))+
#   geom_boxplot()

pdf(paste0(inputPath, "/", aligner, "_align_plots.pdf"), width = 12, height = 12)
par(mfrow=c(3,4))
par(mar=c(2,6,3,1), oma=c(0.5,0.5,3,0.5))
for (variable in 1:12){
  boxplot(x=bam_stats_variables[ ,variable],
          main = title[variable],
          ylab = ylabel[variable],
          col=color[variable],
          xaxt = "n", axes=T, frame.plot = FALSE,
          cex.main=1.3, cex.sub=1.3, cex.lab=1.4, cex.axis=1.3)
}
mtext(bquote(bold("Aligner: ")~bold(.(aligner))), outer = TRUE, cex = 1.5)
dev.off()


#-------- CREATE HTML WITH FORMATTABLE ----------#
# https://cran.r-project.org/web/packages/formattable/vignettes/formattable-data-frame.html
bam_stats_table <- bam_stats
colnames(bam_stats_table) <- c("Sample",	"Depth_Cov",	"Breadth_Cov",
                               "Depth_Only_Cov_Regions",	"Total_Reads",	"Reads_Mapped",
                               "Reads_Unmapped",	"Reads_Paired",	"Bases_Mapped",
                               "Read_Length",	"Base_Quality",	"Mapped_Reads_Perc",	"Unmapped_Reads_Perc")

bam_stats_table[ ,3] <- percent(bam_stats_table[ ,3]/100, digits = 2, format = "f", big.mark = ",")
bam_stats_table[ ,5] <- format(as.numeric(bam_stats_table[ ,5]), big.mark=",")
bam_stats_table[ ,6] <- format(as.numeric(bam_stats_table[ ,6]), big.mark=",")
bam_stats_table[ ,7] <- format(as.numeric(bam_stats_table[ ,7]), big.mark=",")
bam_stats_table[ ,8] <- format(as.numeric(bam_stats_table[ ,8]), big.mark=",")
bam_stats_table[ ,9] <- format(as.numeric(bam_stats_table[ ,9]), big.mark=",")
bam_stats_table[ ,12] <- percent(bam_stats_table[ ,12]/100, digits = 2, format = "f", big.mark = ",")
bam_stats_table[ ,13] <- percent(bam_stats_table[ ,13]/100, digits = 2, format = "f", big.mark = ",")

align_column=c("l", rep("r", 13))

htmlBamStats <- formattable(bam_stats_table, align=align_column, list(
  Depth_Cov = color_tile("lightcyan", "lightcyan4"),
  Breadth_Cov = color_bar("lightblue"),
  Depth_Only_Cov_Regions = color_tile("lightskyblue1", "lightskyblue3"),
  Total_Reads = color_bar("mistyrose"),
  Reads_Mapped = color_bar("lightpink"),
  Reads_Unmapped = color_bar("plum"),
  Reads_Paired = color_bar("khaki"),
  Bases_Mapped = color_bar("navajowhite"),
  Read_Length = color_tile("lightsalmon", "sandybrown"),
  Base_Quality = color_tile("papayawhip", "peachpuff2"),
  Mapped_Reads_Perc = color_bar("lightgreen"),
  Unmapped_Reads_Perc = color_bar("mediumseagreen")
  ))

#---------- SCRIPT TO EXPORT PDF, PNG or JPEG FROM FORMATTABLE ----------#
# https://stackoverflow.com/questions/34983822/how-to-have-r-formattable-rendered-to-pdf-output-and-how-to-have-percents-in-the
# https://github.com/renkun-ken/formattable/issues/26

### It is possible to create exportFormatTable.R file and import the function with source("exportFormatTable.R")
#============================================================================#
# library(htmltools)
# library(webshot)
# suppressMessages(webshot::install_phantomjs())

#' Export a Formattable as PNG, PDF, or JPEG
#'
#' @param f A formattable.
#' @param file Export path with extension .png, .pdf, or .jpeg.
#' @param width Width specification of the html widget being exported.
#' @param height Height specification of the html widget being exported.
#' @param background Background color specification.
#' @param delay Time to wait before taking webshot, in seconds.
#'
#' @importFrom formattable as.htmlwidget
#' @importFrom htmltools html_print
#' @importFrom webshot webshot
#'
#' @export
export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}
#============================================================================#

#export_formattable(htmlBamStats, file = paste0(inputPath, "/", aligner, "_align_summary.pdf"))
#export_formattable(htmlBamStats, file = paste0(inputPath, "/", aligner, "_align_summary.png"))
#export_formattable(htmlBamStats, file = paste0(inputPath, "/", aligner, "_align_summary.jpeg"))


#---------- SCRIPT TO EXPORT HTML FROM FORMATTABLE ----------#
# https://bioinfo.iric.ca/create-a-nice-looking-table-using-r/


html_header="
<head> 
<meta charset=\"utf-8\"> 
<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\"> 
<link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css\">
</head>
<body>
"

suppressWarnings(write(paste(html_header, htmlBamStats, sep=""), file = paste0(inputPath, "/", aligner, "_align_summary.html")))


