##############################
### Orthology normalization ###
##############################
main.dir <- "~/Documents/rhythm_detection_benchmark/random_data"


# install.packages("aggregation")
library("aggregation")

args = commandArgs(trailingOnly=TRUE)
species1 <- args[1]
tissue1 <- args[2]
species2 <- args[3]
tissue2 <- args[4]


if (grepl("mouse", species1)){ species1.name <- "mouse" } else { species1.name <- species1 }
if (grepl("mouse", species2)){ species2.name <- "mouse" } else { species2.name <- species2 }

if (tissue1 != tissue2){
  system("echo You are comparing 2 different tissues")
}


pvalues <- c("raw.pvalue", "default.pvalue", "BH.Q")
for (a in 1:length(pvalues)) {
  p.val <- pvalues[a]
  

  file1.dir <- paste(paste(main.dir, "DATA", species1, tissue1, sep = "/"), "/", sep="")
  file1 <- paste(file1.dir, "normalized_", p.val, "/normalized_pvalue_per_gene.txt", sep = "")
  file2.dir <- paste(paste(main.dir, "DATA", species2, tissue2, sep = "/"), "/", sep="")
  file2 <- paste(file2.dir, "normalized_", p.val, "/normalized_pvalue_per_gene.txt", sep = "")
  
  # check if working directory exist
  if (!file.exists(file1) || !file.exists(file2)){
    system("echo Directory or file does not exist. For Mouse please write: mouse_RNAseq or mouse_microarray. The file expected should be: DATA/species/tissue/normalized_*/normalized_pvalue_per_gene.txt")
    stop()
  }
  
  
  
  # normalization algorithms : PARAMETERS 
  selected.algorithms.for.normalization <- read.csv(paste(main.dir, "scripts/normalization_parameters.csv", sep="/"), nrows = 2, sep = ",") # which type of normalization for genes with several data (microarray ProbIDs or TranscriptIDs)
  per.gene.normalization.algorithm.algorithmname <- as.character(selected.algorithms.for.normalization[selected.algorithms.for.normalization$normalization_step == "per_gene", "normalization_method"])
  per.gene.normalization.algorithm <- paste("pvalue", per.gene.normalization.algorithm.algorithmname, sep="_")
  orthologs.group.normalization.algorithmname <- as.character(selected.algorithms.for.normalization[selected.algorithms.for.normalization$normalization_step == "per_orthologs_group", "normalization_method"])
  if (orthologs.group.normalization.algorithmname == "fisher") {
    orthologs.group.normalization.algorithm <- function(x){ return( fisher(x) ) }
  } else if (orthologs.group.normalization.algorithmname == "sidak") {
    orthologs.group.normalization.algorithm <- function(x){ return( sidak(x) ) }
  } else if (orthologs.group.normalization.algorithmname == "min(bonferroni_pvalues)") {
    orthologs.group.normalization.algorithm <- function(x){ return( min(p.adjust(x, method = "bonferroni", n = length(x))) )}
    orthologs.group.normalization.algorithmname.tmp <- "minimum of bonferroni_pvalues"
  } else { stop("Select a algorithm to normalize P-values of multiple orthologs of species 2 : sidak or fisher or min(bonferroni_pvalues) in scripts/normalization_parameters.csv file") }
        
  system(paste("echo The algorithm selected to normalize p-values per gene is :", per.gene.normalization.algorithm.algorithmname))
  system(paste("echo The algorithm selected to normalize p-values per orthologs group of species2 is :", orthologs.group.normalization.algorithmname.tmp))
  system("echo Modify the scripts/normalization_parameters.csv file if you want to change these p-values or orthologs group normalization algorithms")
  
  
  
  ### Species1
  pvalue.per.gene.1 <- read.table(file1, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  # If initial normalization per gene was needed : we keep only the one selected by pvalue.per.gene.normalization.algorithm argument
  if (length(colnames(pvalue.per.gene.1))>3){
    pvalue.per.gene.1 <- pvalue.per.gene.1[,c("ID", per.gene.normalization.algorithm, "algorithm")]
  }
  colnames(pvalue.per.gene.1) <- c(paste(species1.name, "ID", sep="_"), paste(species1.name, "pvalue", sep="_"), paste(species1.name, "algorithm", sep="_"))
  
  ### Species2
  # If initial normalization per gene was needed : we keep only the one selected by pvalue.per.gene.normalization.algorithm argument
  pvalue.per.gene.2 <- read.table(file2, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if (length(colnames(pvalue.per.gene.2))>3){
    pvalue.per.gene.2 <- pvalue.per.gene.2[,c("ID", per.gene.normalization.algorithm, "algorithm")]
  }
  colnames(pvalue.per.gene.2) <- c(paste(species2.name, "ID", sep="_"), paste(species2.name, "pvalue", sep="_"), paste(species2.name, "algorithm", sep="_"))
  
  
  
  ### ORTHOLOGY ###
  orthologs.file1 <- paste(main.dir, "/DATA/ORTHOLOGY/oma/", species1.name, "_", species2.name, "_oma.txt", sep="")
  orthologs.file2 <- paste(main.dir, "/DATA/ORTHOLOGY/oma/", species2.name, "_", species1.name, "_oma.txt", sep="")
  if (file.exists(orthologs.file1)){
    orthologs.file <- orthologs.file1
  } else if (file.exists(orthologs.file2)){
    orthologs.file <- orthologs.file2
  } else { stop(paste(species1.name, species2.name, "OMA orthology file does not exist.", sep = " ")) }
  
  orthologs <- read.table(orthologs.file, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE, fill=TRUE)
  
  
  ### Then we merge species1 - species2 and their orthology relationships, keeping all informations 
  merged.file <- merge(pvalue.per.gene.1, orthologs, by=paste(species1.name, "ID", sep="_"), all.x = TRUE)
  merged.file <- merge(merged.file, pvalue.per.gene.2, by=paste(species2.name, "ID", sep="_"), all.x = TRUE)
  merged.file <- merged.file[, colnames(merged.file) != paste(species2.name, "ID", sep="_")]
  
  
  ### We separate the merged dataframe in :
  # 1) notconserved.genes = species1 genes with no orthologs in species2
  # 2) conserved.genes.with.uniq.ortholog = Orthologs 1:1 or n:1
  # 3) conserved.genes.with.multiple.orthologs = Orthologs 1:m or n:m
  # 4) conserved.genes.with.multiple.orthologs.with.no.data = Orthologs 1:m or n:m for which species2 orthologs does not have data -> Will have NA for species2 data but keep the orthology info
  notconserved.genes <- unique(subset(merged.file, is.na(merged.file$orthotype)))
  conserved.genes.with.uniq.ortholog <- unique(subset(merged.file, orthotype=="1:1" | orthotype=="n:1"))
  conserved.genes.with.multiple.orthologs <- unique(subset(merged.file, orthotype=="1:m" | orthotype=="n:m"))
  # a unique p-value is affected for the orthologs group "m" using the algorithm selected in scripts/normalization_parameters.csv file.
  conserved.genes.with.multiple.orthologs.normalized <- aggregate(conserved.genes.with.multiple.orthologs[ , paste(species2.name, "pvalue", sep="_")], 
                                                                         list(conserved.genes.with.multiple.orthologs[ , paste(species1.name, "ID", sep="_")], conserved.genes.with.multiple.orthologs[ , paste(species2.name, "algorithm", sep="_")], conserved.genes.with.multiple.orthologs[ , paste(species1.name, "algorithm", sep="_")]),
                                                                         FUN = orthologs.group.normalization.algorithm, drop=FALSE)
  colnames(conserved.genes.with.multiple.orthologs.normalized) <- c(paste(species1.name, "ID", sep="_"), paste(species2.name, "algorithm", sep="_"), paste(species1.name, "algorithm", sep="_"), paste(species2.name, "pvalue", sep="_"))
  
  genes.list.tmp <- unique(conserved.genes.with.multiple.orthologs.normalized[, paste(species1.name, "ID", sep="_")])
  merged.tmp <- unique(merged.file[ merged.file[ , paste(species1.name, "ID", sep="_")] %in% genes.list.tmp, c(paste(species1.name, "ID", sep="_"), "orthotype", "OMAgroup", paste(species1.name, "pvalue", sep="_"), paste(species1.name, "algorithm", sep="_")) ])
  merged.tmp$OMAgroup <- NA
  merged.tmp <- unique(merged.tmp)
  
  conserved.genes.with.multiple.orthologs.normalized <- merge(conserved.genes.with.multiple.orthologs.normalized, merged.tmp, by=c(paste(species1.name, "ID", sep="_"), paste(species1.name, "algorithm", sep="_")), all.y = TRUE)
  conserved.genes.with.multiple.orthologs.with.no.data <- unique(conserved.genes.with.multiple.orthologs[!conserved.genes.with.multiple.orthologs[ , paste(species1.name, "ID", sep="_")] %in% conserved.genes.with.multiple.orthologs.normalized[ , paste(species1.name, "ID", sep="_")],])
  
  full.dataframe <- rbind.data.frame(notconserved.genes, conserved.genes.with.uniq.ortholog, conserved.genes.with.multiple.orthologs.normalized, conserved.genes.with.multiple.orthologs.with.no.data)
  
  # write this final dataframe and keep it as R_object also
  outdir <- paste(main.dir, "/DATA/ORTHOLOGY/result/", sep="")
  outname <- paste(p.val, species1, tissue1, species2, tissue2, sep="_")
  
  if (!dir.exists(outdir)){
    dir.create(file.path(outdir))
  }
  
  write.table(full.dataframe, paste(outdir, outname, ".txt", sep=""), row.names = F, quote = F, sep = "\t")
  
}
