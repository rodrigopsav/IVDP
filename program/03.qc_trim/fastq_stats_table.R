############################################
########## FASTQ STATISTICS TABLE ##########
############################################
my_packages <- c("data.table", "bit64", "dplyr", "formattable", "htmltools", "webshot") # Specify your packages
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])] # Extract not installed packages
if(length(not_installed)) install.packages(not_installed, repos = "https://cran.microsoft.com/", dependencies=TRUE)

suppressMessages(library(data.table))
suppressMessages(library(bit64))
suppressMessages(library(dplyr))
library(formattable)
library(htmltools)
library(webshot)
suppressMessages(webshot::install_phantomjs())

inputPath=inputPath
setwd(inputPath)

#-------- CREATE HTML WITH FORMATTABLE ----------#
# https://cran.r-project.org/web/packages/formattable/vignettes/formattable-data-frame.html
fastq_stats <- fread("fastq_stats.txt", header = T)
fastq_stats <- as.data.frame(fastq_stats)
fastq_stats[,1] <- as.character(fastq_stats[,1])
fastq_stats[,2] <- as.numeric(gsub("\\s*\\([^\\)]+\\)","",as.character( fastq_stats[,2] )))
fastq_stats[,3] <- as.numeric(gsub("\\s*\\([^\\)]+\\)","",as.character( fastq_stats[,3] )))
fastq_stats[,4] <- as.numeric(gsub("[\\%,]", "", fastq_stats[,4]))
fastq_stats[,5] <- as.numeric(gsub("[\\%,]", "", fastq_stats[,5]))
fastq_stats[,6] <- as.numeric(gsub("[\\%,]", "", fastq_stats[,6]))
fastq_stats[,7] <- as.numeric(gsub("[\\%,]", "", fastq_stats[,7]))

### Create format for each column
fastq_stats[ ,2] <- format(as.numeric(fastq_stats[ ,2]), big.mark=",")
fastq_stats[ ,3] <- format(as.numeric(fastq_stats[ ,3]), big.mark=",")
fastq_stats[ ,4] <- percent(fastq_stats[ ,4]/100, digits = 2, format = "f", big.mark = ",")
fastq_stats[ ,5] <- percent(fastq_stats[ ,5]/100, digits = 2, format = "f", big.mark = ",")
fastq_stats[ ,6] <- percent(fastq_stats[ ,6]/100, digits = 2, format = "f", big.mark = ",")
fastq_stats[ ,7] <- percent(fastq_stats[ ,7]/100, digits = 2, format = "f", big.mark = ",")

align_column=c("l", rep("r", 6))

htmlFastqStats <- formattable(fastq_stats, align=align_column, list(
  "Total Reads" = color_bar("mistyrose"),
  "Total Bases" = color_bar("lightpink"),
  "N Bases" = color_bar("lightblue"),
  Q20 = color_bar("lightgreen"),
  Q30 = color_bar("mediumseagreen"),
  GC = color_bar("khaki")
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


#export_formattable(htmlFastqStats, "fastq_summary.pdf")
#export_formattable(htmlFastqStats, "fastq_summary.png")
#export_formattable(htmlFastqStats, "fastq_summary.jpeg")

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

suppressWarnings(write(paste(html_header, htmlFastqStats, sep=""), "fastq_summary.html"))



