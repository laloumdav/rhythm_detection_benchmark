######################################
### Nb of rhythmic genes detected ####
######################################
#install.packages("venn")
library(venn)
library(UpSetR)
library(wesanderson)
library(plyr)
library(maditr)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)
species1 <- args[1]
tissue1 <- args[2]
species2 <- args[3]
tissue2 <- args[4]
species2.algorithm <- args[5]
species2.threshold = args[6]
species2.pvalue <- args[7]
cutoff.genes = as.numeric(args[8])
only.within.conserved.genes <- as.logical(args[9])
only.ortho.one.to.one <- as.logical(args[10])


main.dir <- "~/Documents/rhythm_detection_benchmark"

if (grepl("mouse", species1)){ species1.name <- "mouse" } else { species1.name <- species1 }
if (grepl("mouse", species2)){ species2.name <- "mouse" } else { species2.name <- species2 }


p.val <- "default.pvalue"



### Species1 ### => With species1.pvalue
file.name <- paste(p.val, species1, tissue1, species2, tissue2, sep="_")
orthology.per.gene <- read.table(paste(main.dir, "/DATA/ORTHOLOGY/result/", file.name, ".txt", sep=""), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
### Species2 ### => Retreive species2.pvalue to apply the species2.threshold to these p-values
# Indeed, you can choose to get rhythmic orthologs with raw.pvalues < species2.threshold or default.pvalues < species2.threshold 
file.name.tmp <- paste(species2.pvalue, species1, tissue1, species2, tissue2, sep="_")
orthology.per.gene.tmp <- read.table(paste(main.dir, "/DATA/ORTHOLOGY/result/", file.name.tmp, ".txt", sep=""), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)

orthology.per.gene[[paste(species2.name, "pvalue", sep="_")]] <- orthology.per.gene.tmp[[paste(species2.name, "pvalue", sep="_")]]


# List of algorithms for which there is data :
algorithms.list <- unique(orthology.per.gene[[paste(species1.name, "algorithm", sep="_")]])
algorithms.list <- algorithms.list[!is.na(algorithms.list)]

# Retrieve pvalues of species1 for each algorithm :
species1.pvalues <- dcast(orthology.per.gene,  
                          orthology.per.gene[[paste(species1.name, "ID", sep="_")]] ~ orthology.per.gene[[paste(species1.name, "algorithm", sep="_")]], 
                          value.var = paste(species1.name, "pvalue", sep="_"), 
                          fun.aggregate = mean, drop=TRUE, fill=1)
colnames(species1.pvalues)[1] <- paste(species1.name, "ID", sep="_")
# If not previously removed, we replace NA by value=1 because it correspond to outsider genes with absolutly NO signal of rhythmicity
species1.pvalues <- replace(species1.pvalues, is.na(species1.pvalues), 1)




# Conservation info:
conservation <- unique(orthology.per.gene[ , c(paste(species1.name, "ID", sep="_"), "orthotype")])
#### Orthology ####
orthologs.pvalues <- dcast(orthology.per.gene,  
                           orthology.per.gene[, paste(species1.name, "ID", sep="_")] ~ orthology.per.gene[, paste(species2.name, "algorithm", sep="_")],
                           value.var = paste(species2.name, "pvalue", sep="_"), 
                           fun.aggregate = mean, drop=TRUE, fill=NULL)
# Remove NA column if exist
orthologs.pvalues <- orthologs.pvalues[, colnames(orthologs.pvalues)!="NA"]
column.names <- colnames(orthologs.pvalues)
column.names <- c(paste(species1.name, "ID", sep="_"), paste(species2.name, column.names[-1], sep="_"))
colnames(orthologs.pvalues) <- column.names
# Is there data for the ortholog gene in species2 ?
orthologs.pvalues$species2_data_available <- apply(orthologs.pvalues[,-1], 1, function(x) all(!is.nan(x)))
# Conserved genes: Affect TRUE if the gene of species1 has one or more orthologs in species2
orthologs.pvalues <- merge(orthologs.pvalues, conservation, by=paste(species1.name, "ID", sep="_"))
orthologs.pvalues$conserved_gene <- !is.na(orthologs.pvalues$orthotype)


# Remove conserved genes that do not have data in species2 :
master <- subset(orthologs.pvalues, !(species2_data_available==FALSE & conserved_gene==TRUE))
# retrieve species1 data
master <- merge(master, species1.pvalues, by=paste(species1.name, "ID", sep="_"))


###
### Tot number of orthologs 
###
sub.master <- master
# If argument only.ortho.one.to.one = TRUE : we keep only orthologs 1:1
if (isTRUE(only.ortho.one.to.one)){ sub.master <- subset(sub.master, orthotype == "1:1") }

# For each algorithm we keep the first "cutoff.genes" orthologs (ranked by p-values) in species 1
# And keep those which are rhythmic in species 2 as well:
list.rhythmic.ortho <- list()
for (j in 1:length(algorithms.list)) {
  data <- sub.master[order(sub.master[, algorithms.list[j]]), ]
  # Keep cutoff.genes first genes:
  data <- data[1:cutoff.genes, ]
  # Keep genes which are rhythmic in species 2 as well:
  data <- subset(data, data[[paste(species2.name, species2.algorithm, sep="_")]] <= species2.threshold)
  
  list.rhythmic.ortho[[algorithms.list[j]]] <- unique(data[[paste(species1.name, "ID", sep="_")]])
}

# Take all genes considered
genes <- unique(unlist(list.rhythmic.ortho))

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

intersect.table.rhythmic.ortho <- FromList(list.rhythmic.ortho, genes)

#intersect.table.rhythmic.ortho <- intersect.table.rhythmic.ortho[, c("ARS", "GeneCycle", "empJTK", "elements")]


################
#  Upset plot  #
################
colors.list <- c(wes_palette("Darjeeling1"), wes_palette("GrandBudapest2"))

# function to draw colored histogram corresponding with the smaller p-value threshold
upsetFunc <- function(input.data) {
  data <- (input.data["comb.retrieved.smaller.threshold"] == "yes")
}
# metadata to get colored line for each algorithm
metadata <- data.frame(sets=names(list.rhythmic.ortho), algorithm=names(list.rhythmic.ortho))

intersect.table.rhythmic.ortho <- intersect.table.rhythmic.ortho[, colnames(intersect.table.rhythmic.ortho) %in% c("LS", "JTK", "empJTK", "RAIN", "GeneCycle", "meta2d", "ARS")]
intersect.table.rhythmic.ortho$algo.nb <- apply(intersect.table.rhythmic.ortho, 1, sum)

upset.query <- list(query = intersects, 
                      params = list("LS", "JTK", "empJTK", "RAIN", "GeneCycle", "meta2d", "ARS"), 
                      color = "cyan4", active = F, 
                      query.name = "intersection with all algorithms")

pfd.name <- paste("upset_RhythmicOrtho_", p.val, "_", species1.name, "_", species2.name, "_", tissue1, ".pdf", sep="")
pfd.name <- paste(main.dir, "/RESULTS/proportion_rhythmic_orthologs_images/", pfd.name, sep="")
pdf(file = pfd.name, onefile=FALSE) # or other device
print(
  upset(intersect.table.rhythmic.ortho, 
      query.legend = "top",
      mainbar.y.label = paste("Intersection size", sep=""),
      nsets=length(list.rhythmic.ortho), 
      nintersects=60,
      order.by = c("freq", "degree"), decreasing = c(TRUE,TRUE),
      set.metadata = list(data = metadata, plots = list(list(type = "matrix_rows", column = "algorithm", colors = c(JTK=colors.list[1], LS=colors.list[6], empJTK=colors.list[7], ARS=colors.list[4], RAIN=colors.list[5], meta2d=colors.list[2], GeneCycle="blue"), alpha = 0.5))),
      queries = list(upset.query))
)
dev.off()


