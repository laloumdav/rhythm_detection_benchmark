######################################
### Nb of rhythmic genes detected ####
######################################
#install.packages("venn")
library(venn)
library(UpSetR)
library(wesanderson)
library(plyr)

main.dir <- "~/Documents/rhythm_detection_benchmark"

args = commandArgs(trailingOnly=TRUE)
species <- args[1]
tissue <- args[2]
threshold <- as.numeric(args[3])


if (grepl("mouse", species)){ species.name <- "mouse" } else { species.name <- species }

# normalization methods : PARAMETERS 
selected.methods.for.normalization <- read.csv(paste(main.dir, "scripts/normalization_parameters.csv", sep="/"), nrows = 2, sep = ",") # which type of normalization for genes with several data (microarray ProbIDs or TranscriptIDs)
per.gene.normalization.method.name <- as.character(selected.methods.for.normalization[selected.methods.for.normalization$normalization_step == "per_gene", "normalization_method"])
per.gene.normalization.method <- paste("pvalue", per.gene.normalization.method.name, sep="_")


file.dir <- paste(paste(main.dir, "DATA", species, tissue, sep = "/"), "/", sep="")

pvalues <- c("raw.pvalue", "default.pvalue", "BH.Q")
for (a in 1:length(pvalues)) {
  p.val <- pvalues[a]
  file.name <- paste(file.dir, "normalized_", p.val, "/normalized_pvalue_per_gene.txt", sep = "")
  
  normalized.file <- read.table(file.name, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  # If initial normalization per gene was needed : we keep only the one selected by pvalue.per.gene.normalization.method argument
  if (length(colnames(normalized.file))>3){
    normalized.file <- normalized.file[,c("ID", per.gene.normalization.method, "algorithm")]
  }
  colnames(normalized.file) <- c("ID", "pvalue", "algorithm")
  
  
  tot.genes.nb = length(unique(normalized.file$ID))
  list.genes.per.algorithm <- split(normalized.file, f = normalized.file$algorithm)
  
  # keep first rhythmic genes list : first "threshold" genes.
  list.rhythmicgenes.per.algorithm_threshold <- list()
  for (i in 1:length(list.genes.per.algorithm)){
    data <- list.genes.per.algorithm[[i]]
    
    algorithm.name <- data[1, "algorithm"]
    # subset of first "threshold" rhythmic genes
    data <- data[order(data[, "pvalue"]), ]
    data <- data[1:threshold, ]
    list.rhythmicgenes.per.algorithm_threshold[[algorithm.name]] <- unique(data$ID)
  }
  
  # Take genes kept with the highest threshold
  genes <- unique(unlist(list.rhythmicgenes.per.algorithm_threshold))
  
  # Function to unlist intersections of rhythmic genes between algorithms which keep genes ID 
  FromList <- function(input.list, elements) {
    data <- unlist(lapply(input.list, function(x) { x <- as.vector(match(elements, x)) } ))
    data[is.na(data)] <- as.integer(0)
    data[data != 0] <- as.integer(1)
    data <- data.frame(matrix(data, ncol = length(input.list), byrow = F))
    names(data) <- names(input.list)
    data$elements <- elements
    return(data)
  }
  
  intersect.table.rhythmicgenes_threshold <- FromList(list.rhythmicgenes.per.algorithm_threshold, genes)
  intersect.table.rhythmicgenes <- intersect.table.rhythmicgenes_threshold
  
  
  #intersect.table.rhythmicgenes <- intersect.table.rhythmicgenes[, c("ARS", "GeneCycle", "empJTK", "elements")]
  
  
  ################
  #  Upset plot  #
  ################
  colors.list <- c(wes_palette("Darjeeling1"), wes_palette("GrandBudapest2"))
  
  # function to draw colored histogram corresponding with the smaller p-value threshold
  upsetFunc <- function(input.data) {
    data <- (input.data["comb.retrieved.smaller.threshold"] == "yes")
  }
  # metadata to get colored line for each algorithm
  metadata <- data.frame(sets=names(list.genes.per.algorithm), algorithm=names(list.genes.per.algorithm))
  
  simple.intersect.table.rhythmicgenes <- intersect.table.rhythmicgenes[, colnames(intersect.table.rhythmicgenes) %in% c("LS", "JTK", "empJTK", "RAIN", "GeneCycle", "meta2d", "ARS")]
  simple.intersect.table.rhythmicgenes$algo.nb <- apply(simple.intersect.table.rhythmicgenes, 1, sum)
  count.algo.nb <- plyr::count(simple.intersect.table.rhythmicgenes$algo.nb)
  
  
  pfd.name <- paste("orderedGenes_upset_DEGREE_", p.val, "_", species, "_", tissue, ".pdf", sep="")
  pfd.name <- paste(main.dir, "/RESULTS/UpSet_images/orderedGenes_UpSet_images/", pfd.name, sep="")
  pdf(file = pfd.name, onefile=FALSE) # or other device
  print(
    upset(intersect.table.rhythmicgenes, 
        query.legend = "top",
        mainbar.y.label = paste("Intersection size \n cutoff=", threshold, sep=""),
        nsets=length(list.genes.per.algorithm), 
        nintersects=60,
        order.by = c("freq", "degree"), decreasing = c(TRUE,TRUE),
        set.metadata = list(data = metadata, plots = list(list(type = "matrix_rows", column = "algorithm", colors = c(JTK=colors.list[1], LS=colors.list[6], empJTK=colors.list[7], ARS=colors.list[4], RAIN=colors.list[5], meta2d=colors.list[2], GeneCycle="blue"), alpha = 0.5))))
  )
  dev.off()
  
  
}


