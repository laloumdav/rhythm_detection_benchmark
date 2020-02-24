#####################
#### GeneCycle ######
#####################

main.dir <- "~/Documents/rhythm_detection_benchmark/random_data"

install.packages(paste(main.dir, "scripts/GeneCycle", sep = "/"), repos = NULL, type = "source")
library(GeneCycle)
library(readr)
library(plyr)


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
if ( length(duplicated(time.points) == FALSE) != length(duplicated(time.points) == FALSE) ){
  print("ERROR: The number of replicates isn't identical between the timepoints")
} else { replicates.number <- as.numeric(unique(count(time.points)$freq)) }


# If there is several replicates per time-point, we are going to transform replicates as new days measurements 
# Indeed, robust.spectrum function of GeneCycle can not deal with several replicates per time-point

if ( any(duplicated(time.points)) == TRUE ) {
  GeneCycle.time.points <- time.points[duplicated(time.points)==FALSE]
  nb.of.replicates <- unique(plyr::count(time.points)$freq)
  GeneCycle.seq <- seq(from=2, to=ncol(raw.dataset), by=nb.of.replicates)
  # interval.time
  interval.time <- GeneCycle.time.points[3]-GeneCycle.time.points[2]
  
  GeneCycle.raw.dataset <- raw.dataset[ , c(1, GeneCycle.seq)]
  for (i in 1:(nb.of.replicates-1)) {
    GeneCycle.raw.dataset <- cbind(GeneCycle.raw.dataset, raw.dataset[ , c(GeneCycle.seq+i)])
  }
  
  unique.time.points <- seq(from = min(GeneCycle.time.points), by = interval.time, length.out = length(time.points))
  length(unique.time.points)
  
  ### GeneCycle :
  raw.dataset.IDs <- as.data.frame(GeneCycle.raw.dataset$ID)
  raw.dataset.transposed <- t(GeneCycle.raw.dataset[,-1])
  
} else {
  # interval.time
  unique.time.points <- unique(time.points)
  
  ### GeneCycle :
  raw.dataset.IDs <- as.data.frame(raw.dataset$ID)
  raw.dataset.transposed <- t(raw.dataset[,-1])
}

raw.dataset.transposed <- raw.dataset.transposed+0.0001
## Lets use the robust regression based approach (Ahdesmaki et al. 2007) 
##with robust.spectrum function. Because we known the periodicity time.
## (computation will take a lot of time depending on how many permutations are used per time series and time series length).
GeneCycle.output <- robust.spectrum(x = raw.dataset.transposed,
                                    t = unique.time.points,
                                    periodicity.time = 24,
                                    #noOfPermutations = 50,
                                    algorithm = "regression")

GeneCycle <- data.frame(ID = raw.dataset.IDs, 
                        raw.pvalue = GeneCycle.output, 
                        default.pvalue = GeneCycle.output, 
                        BH.Q = p.adjust(GeneCycle.output, method = "BH"), 
                        period = NA, phase = NA, amplitude = NA)

colnames(GeneCycle)[1] <- "ID"


write.table(GeneCycle, paste(file.dir, "GeneCycle.txt", sep = ""), row.names = F, quote = F, sep = "\t")

