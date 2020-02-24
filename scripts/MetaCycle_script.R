#################
### METACYCLE ###
#################
main.dir <- "~/Documents/rhythm_detection_benchmark"

library(readr)
library(plyr)
## Load modified MetaCycle Package :
install.packages(paste(main.dir, "scripts/MetaCycle", sep = "/"), repos = NULL, type = "source")
library(MetaCycle)


args = commandArgs(trailingOnly=TRUE)
species <- args[1]
tissue <- args[2]
tissue.file <- paste(tissue, ".txt", sep="")


file.dir <- paste(paste(main.dir, "DATA", species, tissue, sep = "/"), "/", sep="")

# check if working directory exist
if (dir.exists(file.dir)==FALSE){
  stop("Directory/file doesn't exist. Maybe: mouse_RNAseq or mouse_microarray ???")
}


#### SCRIPT ####
system(paste("echo Running the file :", tissue.file, sep=" "))

file.name <- paste(file.dir, tissue.file, sep = "")
raw.dataset <- read.table(file.name, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)

time.points <- colnames(raw.dataset)[-1]
time.points <- parse_number(time.points)

# If there is several replicates per time-point, we are going to keep only one of them for running ARS algorithm
# Indeed, ARS can not deal with several replicates per time-point

if ( any(duplicated(time.points)) == TRUE ) {
  ARS.time.points <- time.points[duplicated(time.points)==FALSE]
  nb.of.replicates <- unique(plyr::count(time.points)$freq)
  ARS.seq <- seq(from=2, to=ncol(raw.dataset), by=nb.of.replicates)
  
  ARS.raw.dataset <- raw.dataset[ , c(1, ARS.seq)]
  write.table(ARS.raw.dataset, paste(file.dir, "ARS.raw.dataset.txt", sep = ""), row.names = F, quote = F, sep = "\t")
}


### MetaCycle :
if ( any(duplicated(time.points)) == TRUE ) {
  ARS.output.default <- meta2d(infile=paste(file.dir, "ARS.raw.dataset.txt", sep = ""), filestyle="txt",
                                    outputFile=FALSE,
                                    minper = 20, maxper = 28,
                                    timepoints=ARS.time.points,
                                    cycMethod="ARS", outRawData=FALSE, 
                                    raw.pVal=FALSE) 
  ARS.output.raw <- meta2d(infile=paste(file.dir, "ARS.raw.dataset.txt", sep = ""), filestyle="txt",
                                 outputFile=FALSE,
                                 minper = 20, maxper = 28,
                                 timepoints=ARS.time.points,
                                 cycMethod="ARS", outRawData=FALSE, 
                                 raw.pVal=TRUE) ### NEW OPTION ADDED TO GET RAW P-VALUES AS OUTPUT ###
  
  ARS.default <- ARS.output.default$meta
  ARS.raw <- ARS.output.raw$meta
}

metaCycle.output.default <- meta2d(infile=file.name, filestyle="txt",
                             outputFile=FALSE,
                             minper = 20, maxper = 28,
                             timepoints=time.points,
                             cycMethod=c("LS", "JTK", "ARS"), outRawData=FALSE, 
                             raw.pVal=FALSE) ### NEW OPTION ADDED TO GET RAW P-VALUES AS OUTPUT ###
metaCycle.output.raw <- meta2d(infile=file.name, filestyle="txt",
                             outputFile=FALSE,
                             minper = 20, maxper = 28,
                             timepoints=time.points,
                             cycMethod=c("LS", "JTK", "ARS"), outRawData=FALSE, 
                             raw.pVal=TRUE) ### NEW OPTION ADDED TO GET RAW P-VALUES AS OUTPUT ###


meta.raw <- metaCycle.output.raw$meta
meta.default <- metaCycle.output.default$meta

if ( any(duplicated(time.points)) == TRUE ) {
  meta.default <- cbind(meta.default, ARS.default[, grep("ARS", colnames(ARS.default))])
  meta.raw <- cbind(meta.raw, ARS.raw[, grep("ARS", colnames(ARS.raw))])
  
  # remove the created file specialy for ARS : "ARS.raw.dataset.txt"
  file.remove( paste(file.dir, "ARS.raw.dataset.txt", sep = "") )
}  

if (grepl("ARS", do.call("paste",as.list(colnames(meta.raw))))){
ARS <- meta.raw[ , c("CycID", "ARS_pvalue")]
ARS <- cbind(ARS, meta.default[ , c("ARS_pvalue", "ARS_period", "ARS_adjphase", "ARS_amplitude")])
colnames(ARS) <- c("ID", "raw.pvalue", "default.pvalue", "period", "phase", "amplitude")
}
if (grepl("JTK", do.call("paste",as.list(colnames(meta.raw))))){
JTK <- meta.raw[ , c("CycID", "JTK_pvalue")]
JTK <- cbind(JTK, meta.default[ , c("JTK_pvalue", "JTK_period", "JTK_adjphase", "JTK_amplitude")])
colnames(JTK) <- c("ID", "raw.pvalue", "default.pvalue", "period", "phase", "amplitude")
}
if (grepl("LS", do.call("paste",as.list(colnames(meta.raw))))){
LS <- meta.raw[ , c("CycID", "LS_pvalue")]
LS <- cbind(LS, meta.default[ , c("LS_pvalue", "LS_period", "LS_adjphase", "LS_amplitude")])
colnames(LS) <- c("ID", "raw.pvalue", "default.pvalue", "period", "phase", "amplitude")
}
if (grepl("meta2d", do.call("paste",as.list(colnames(meta.raw))))){
meta2d <- meta.raw[ , c("CycID", "meta2d_pvalue")]
meta2d <- cbind(meta2d, meta.default[ , c("meta2d_pvalue", "meta2d_period", "meta2d_phase", "meta2d_AMP")])
colnames(meta2d) <- c("ID", "raw.pvalue", "default.pvalue", "period", "phase", "amplitude")
}




if (exists("ARS")){ write.table(ARS, paste(file.dir, "ARS.txt", sep = ""), row.names = F, quote = F, sep = "\t") }
if (exists("JTK")){ write.table(JTK, paste(file.dir, "JTK.txt", sep = ""), row.names = F, quote = F, sep = "\t") }
if (exists("LS")){ write.table(LS, paste(file.dir, "LS.txt", sep = ""), row.names = F, quote = F, sep = "\t") }
if (exists("meta2d")){ write.table(meta2d, paste(file.dir, "meta2d.txt", sep = ""), row.names = F, quote = F, sep = "\t") }


