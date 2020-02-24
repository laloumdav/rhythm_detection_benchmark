#################
#### RAIN #######
#################
main.dir <- "~/Documents/rhythm_detection_benchmark"

library(readr)
library(plyr)
## Load modified rain Package : 
install.packages(paste(main.dir, "scripts/rain", sep = "/"), 
                 repos = NULL, 
                 type = "source")
library(rain)


args = commandArgs(trailingOnly=TRUE)
species <- args[1]
tissue <- args[2]


tissue.file <- paste(tissue, ".txt", sep="")
file.dir <- paste(paste(main.dir, "DATA", species, tissue, sep = "/"), "/", sep="")

# check if working directory exist
if (dir.exists(file.dir)==FALSE){
  stop("Directory/file doesn't exist. Maybe: mouse_RNAseq or mouse_microarray ???")
}


### SCRIPT ###
system(paste("echo Running the file :", tissue.file, sep=" "))

file.name <- paste(file.dir, tissue.file, sep = "")
raw.dataset <- read.table(file.name, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)


# time.points
time.points <- colnames(raw.dataset)[-1]
time.points <- parse_number(time.points)
# Number of replicates for each timepoint
if (isUnique(unique(plyr::count(time.points)$freq))==FALSE){
  print("ERROR: The number of replicates isn't identical between the timepoints")
} else { replicates.number <- as.numeric(unique(plyr::count(time.points)$freq))}
# interval.time
unique.time.points <- unique(time.points)
interval.time <- unique.time.points[3]-unique.time.points[2]


### RAIN :
raw.dataset.IDs <- as.data.frame(raw.dataset$ID)
raw.dataset.ready.for.rain <- raw.dataset[,-1]
RAIN.output.default <- rain(t(raw.dataset.ready.for.rain),
                            deltat = interval.time,
                            period = 24, period.delta=4,
                            nr.series = replicates.number,
                            raw.pVal = FALSE,  ### Added option to get raw p-values
                            verbose=TRUE, method='independent')
RAIN.output.raw <- rain(t(raw.dataset.ready.for.rain),
                            deltat = interval.time,
                            period = 24, period.delta=4,
                            nr.series = replicates.number,
                            raw.pVal = TRUE,  ### Added option to get raw p-values
                            verbose=TRUE, method='independent')
RAIN.output.raw <- as.data.frame(RAIN.output.raw$pVal)

RAIN <- cbind(raw.dataset.IDs, RAIN.output.raw, RAIN.output.default, NA)
RAIN <- RAIN[,c(1,2,3,6,4,7)]
colnames(RAIN) <- c("ID", "raw.pvalue", "default.pvalue", "period", "phase", "amplitude")



write.table(RAIN, paste(file.dir, "RAIN.txt", sep = ""), row.names = F, quote = F, sep = "\t")

