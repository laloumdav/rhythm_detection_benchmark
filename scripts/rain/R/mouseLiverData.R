
#' @name menetRNASeqMouseLiver
#' @title Time courses of gene expression in mouse liver
#' @description Temporal gene expression profiling in mouse liver, measured by 
#' high throughput sequencing. Profiles shows the changes in gene expression 
#' in mice under 12 h light: 12 h dark conditions. Data were originally 
#' published by Menet et. 
#' al. (2012) under Creative Commons License 3.0: 
#' \url{http://creativecommons.org/licenses/by/3.0/}. The data are available 
#' in the public domain at GEO \url{http://www.ncbi.nlm.nih.gov/geo/} as 
#' Dataset entry GSE36916.
#' @docType data
#' @usage menetRNASeqMouseLiver
#' @format a \code{data.frame} containing the normalized expression values for 
#' each gene and time point in two repeats. First number shows the time of 
#' measurement in 'Zeitgeber Time' (ZT) whereas ZT_0 is the time of light on. 
#' The second number specifies the biological replicate.
#' 
#' @references Menet, J. S., Rodriguez, J., Abruzzi, K. C., & Rosbash, M. 
#' (2012). Nascent-Seq reveals novel features of mouse circadian 
#' transcriptional regulation. \emph{eLife}, \bold{1(0)}, e00011. 
#' doi:10.7554/eLife.00011
#' 
#' @author Paul F. Thaben \email{paul.thaben@@charite.de}
NULL 
