#######################
### DATA transform ####
#### into GeneID ######
#######################

args = commandArgs(trailingOnly=TRUE)
species <- args[1]
tissue <- args[2]

main.dir <- "~/Documents/rhythm_detection_benchmark"

file.dir <- paste(paste(main.dir, "DATA", species, tissue, sep = "/"), "/", sep="")
tissue.file <- list.files(file.dir, pattern = "raw.txt")

# check if working directory exist
if (dir.exists(file.dir)==FALSE){
  stop("Directory/file doesn't exist. Maybe: mouse_RNAseq or mouse_microarray ???")
}

### Check if unique 
if (length(tissue.file)!=1){
  stop("raw file is not unique or does not exist. Must be named as tissue_raw.txt in DATA/species/tissue/ directory.")
}


IDs.crossRef.dir <- paste(paste(main.dir, "DATA", species, sep = "/"), "/", sep="")


### SCRIPT ###
file.name <- paste(file.dir, tissue.file, sep = "")
raw.dataset <- read.table(file.name, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# If needed, check if IDs.crossRef.file exist in the good directory place and is unique
if (!grepl("ENSDARG|ENSRNOG|ENSMUSG|ENSPANG|AGA|AAE|FBg", raw.dataset[2, 1])){
  if (length(list.files(IDs.crossRef.dir, pattern = "ID")) !=1){
    stop("ERROR: IDs.crossReference_file doesn't exist or isn't unique. Should be a unique file in rhythm_detection_benchmark/DATA/species/  directory")
  }
  
  IDs.crossRef.file <- paste(IDs.crossRef.dir, list.files(IDs.crossRef.dir, pattern = "ID"), sep="")
  IDs.crossRef <- read.table(IDs.crossRef.file, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  
  # To select the good column to be merged (AFFY.ID or Transcript.ID for instance)
  # Check if there is a unique column un IDs.crossRef matching with some IDs of raw data (here, the number displayed by unlist(apply()) is the column number)
  if (length(unique(unlist(apply(IDs.crossRef, 1, function(x) which(x == raw.dataset$ID[1] | x == raw.dataset$ID[3] | x == raw.dataset$ID[5] ))))) != 1){
    print("ERROR: No columns in IDs.crossReference_file matches with IDs of raw data")
    stop()
  }
  # To select the good column name of IDs.crossRef to be merged (the number displayed by unlist(apply()) is the column number)
  column.name <- colnames(IDs.crossRef)[unique(unlist(apply(IDs.crossRef, 1, function(x) which(x == raw.dataset$ID[1] | x == raw.dataset$ID[3] | x == raw.dataset$ID[5] ))))]
  
  # Merging:
  raw.dataset.merged <- merge(IDs.crossRef, raw.dataset, by.x=column.name, by.y="ID", all.x=T)
  
  # Keeping protein_coding_genes:
  if ("Gene.Type" %in% colnames(raw.dataset.merged)){
    raw.dataset.merged <- subset(raw.dataset.merged, Gene.Type=="protein_coding")
  }
  
  # Check for Gene.ID or GeneID entries
  if ("GeneID" %in% colnames(raw.dataset.merged)){
    colnames(raw.dataset.merged)[grep("GeneID", colnames(raw.dataset.merged))] <- "ID"
  }
  
  # THEN: We remove prob.ID if it is affected to several Gene.ID (long process) 
  rows.to.remove <- c()
  for (i in 1:nrow(raw.dataset.merged)) {
    if ((i %in% rows.to.remove)==FALSE) {
      prob.ID <- raw.dataset.merged[i, column.name]
      if (length(unique(raw.dataset.merged$Gene.ID[grep(prob.ID, raw.dataset.merged[, column.name])])) != 1){
        rows.to.remove <- c(rows.to.remove, grep(prob.ID, raw.dataset.merged[, column.name]))
      }
    }
  }
    
  rows.to.remove <- unique(rows.to.remove)
  
  # So, we remove prob.IDs refered to several Gene.IDs:
  if (is.null(rows.to.remove)==FALSE){
    raw.dataset.merged <- raw.dataset.merged[-rows.to.remove, ]
  }
  
  # Then keep only Gene.ID and Data
  colnumbers.to.keep <- grep("Gene.ID|ZT|CT", colnames(raw.dataset.merged))
  raw.dataset.genes <- raw.dataset.merged[ , colnumbers.to.keep]
  colnames(raw.dataset.genes)[1] <- "ID"
  
} else { raw.dataset.genes <- raw.dataset }

# remove .x .y in colnames if there is several replicates: 
colnames(raw.dataset.genes) <- gsub("\\..*", "", colnames(raw.dataset.genes))

# We remove genes with NO expression (=0) at all time-points :
raw.dataset.genes <- raw.dataset.genes[!apply(raw.dataset.genes[,-1], 1, function(x) all(x==0)), ]


# Write final dataset :
write.table(raw.dataset.genes, gsub("_raw", "", file.name), row.names = F, quote = F, sep = "\t")

