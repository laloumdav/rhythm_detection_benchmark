########################################
### Density distribution of pvalues ####
### in different sub-group #############
########################################

args = commandArgs(trailingOnly=TRUE)
species <- args[1]
tissue <- args[2]


main.dir <- "~/Documents/rhythm_detection_benchmark"
main.dir.ramdom.data <- "~/Documents/rhythm_detection_benchmark/random_data"


if (grepl("mouse", species)){ species.name <- "mouse" } else { species.name <- species }

# normalization methods : PARAMETERS 
selected.methods.for.normalization <- read.csv(paste(main.dir, "scripts/normalization_parameters.csv", sep="/"), nrows = 2, sep = ",") # which type of normalization for genes with several data (microarray ProbIDs or TranscriptIDs)
per.gene.normalization.method.name <- as.character(selected.methods.for.normalization[selected.methods.for.normalization$normalization_step == "per_gene", "normalization_method"])
per.gene.normalization.method <- paste("pvalue", per.gene.normalization.method.name, sep="_")
orthologs.group.normalization.method.name <- as.character(selected.methods.for.normalization[selected.methods.for.normalization$normalization_step == "per_orthologs_group", "normalization_method"])
if (orthologs.group.normalization.method.name == "fisher") {
  orthologs.group.normalization.method <- function(x){ return( fisher(x) ) }
} else if (orthologs.group.normalization.method.name == "sidak") {
  orthologs.group.normalization.method <- function(x){ return( sidak(x) ) }
} else if (orthologs.group.normalization.method.name == "min(bonferroni_pvalues)") {
  orthologs.group.normalization.method <- function(x){ return( min(p.adjust(x, method = "bonferroni", n = length(x))) )}
  orthologs.group.normalization.method.name.tmp <- "minimum of bonferroni_pvalues"
} else { stop("Select a algorithm to normalize P-values of multiple orthologs of species 2 : sidak or fisher or min(bonferroni_pvalues) in scripts/normalization_parameters.csv file") }

system(paste("echo The method selected to normalize p-values per gene is :", per.gene.normalization.method))
system(paste("echo The method selected to normalize p-values per orthologs group of species2 is :", orthologs.group.normalization.method.name.tmp))
system("echo Modify the scripts/normalization_parameters.csv file if you want to change these p-values or orthologs group normalization methods")


p.val <- "default.pvalue"

###############################################################
###############################################################
########  1) On original data     #############################
###############################################################
###############################################################

####################
# Global distribution : 
####################
file.dir <- paste(paste(main.dir, "DATA", species, tissue, sep = "/"), "/", sep="")

file <- paste(file.dir, "normalized_", p.val, "/normalized_pvalue_per_gene.txt", sep = "")
normalized.file <- read.table(file, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
# If initial normalization per gene was needed : we keep only the one selected by pvalue.per.gene.normalization.method argument
if (length(colnames(normalized.file))>3){
  normalized.file <- normalized.file[,c("ID", per.gene.normalization.method, "algorithm")]
}
colnames(normalized.file) <- c("ID", "pvalue", "algorithm")
normalized.file$quantile <- rep("full dataset", nrow(normalized.file))

####################
# Through quartiles : 
####################
########################################################################
# 1) the 1st quartile of mediane gene expression level 
# 2) the 4th quartile of mediane gene expression level 
########################################################################
raw.species.file.name <- paste(file.dir, tissue, ".txt", sep = "")
raw.species.data <- read.table(raw.species.file.name, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
raw.species.data$mediane.expr <- apply(raw.species.data[, -1], 1, median)
raw.species.data <- raw.species.data[, c("ID", "mediane.expr")]
# if several transcripts per gene, we keep the mean of mediane gene expr to get 1 value per gene
raw.species.data <- aggregate(raw.species.data["mediane.expr"], raw.species.data[1], mean)

normalized.file <- merge(normalized.file[, -4], raw.species.data, by = "ID")

algo.list <- unique(normalized.file$algorithm)
t.test.table <- data.frame(method=NA, n=NA, Mean)
for (i in 1:length(algo.list)) {
  subset.algo <- subset(normalized.file, algorithm == algo.list[i])
  rhythmic.genes <- subset(subset.algo, pvalue <= 0.005)
  non.rhythmic.genes <- subset(subset.algo, pvalue > 0.01)
  non.rhythmic.genes <- non.rhythmic.genes[sample(1:nrow(non.rhythmic.genes), size = nrow(rhythmic.genes)), ]

  nrow(non.rhythmic.genes)
  t.test(rhythmic.genes$mediane.expr, non.rhythmic.genes$mediane.expr, alternative = "greater")
  
}






