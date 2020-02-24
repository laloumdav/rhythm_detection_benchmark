########################################
### Density distribution of pvalues ####
### in different sub-group #############
########################################
require(RColorBrewer)
library(tidyverse)
library(wesanderson)
library(magick)
library(ggridges)
library(cowplot)



species <- "mouse_microarray"
tissue <- "liver"
limited.time.points.1 <- "exploration_within_restricted_timepoints_zebrafish"
limited.time.points.2 <- "exploration_within_restricted_timepoints"
limited.time.points.3 <- "exploration_within_restricted_timepoints_baboon"


main.dir <- "~/Documents/rhythm_detection_benchmark"



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
########################################
# limited.time.points.0 ####################
########################################
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
### We separate genes en 4 quartiles
raw.species.data <- raw.species.data %>% dplyr::mutate(quantile = ntile(mediane.expr, 4))

normalized.file.quartiles <- merge(normalized.file[, -4], raw.species.data, by = "ID")
# Keep only 1st and 4th quartiles :
normalized.file.quartiles <- subset(normalized.file.quartiles, quantile %in% c(1,4))
normalized.file.quartiles <- normalized.file.quartiles[, c("ID", "pvalue", "algorithm", "quantile")]
normalized.file.quartiles[normalized.file.quartiles$quantile == 4, "quantile"] <- "4th quartile of gene expression full dataset"
normalized.file.quartiles[normalized.file.quartiles$quantile == 1, "quantile"] <- "1st quartile of gene expression full dataset"

normalized.file <- rbind(normalized.file, normalized.file.quartiles)

# To conserve to same order of groups in the plot : 
ordered.group <- c("random data", "4th quartile of gene expression random data", "1st quartile of gene expression random data", "full dataset", "4th quartile of gene expression full dataset", "1st quartile of gene expression full dataset", "known cycling genes")
normalized.file$quantile <- factor(normalized.file$quantile, levels=ordered.group)

normalized.file.0 <- normalized.file

########################################
# limited.time.points.1 ####################
########################################
file.dir <- paste(paste(main.dir, "DATA", species, tissue, limited.time.points.1, sep = "/"), "/", sep="")
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
### We separate genes en 4 quartiles
raw.species.data <- raw.species.data %>% dplyr::mutate(quantile = ntile(mediane.expr, 4))

normalized.file.quartiles <- merge(normalized.file[, -4], raw.species.data, by = "ID")
# Keep only 1st and 4th quartiles :
normalized.file.quartiles <- subset(normalized.file.quartiles, quantile %in% c(1,4))
normalized.file.quartiles <- normalized.file.quartiles[, c("ID", "pvalue", "algorithm", "quantile")]
normalized.file.quartiles[normalized.file.quartiles$quantile == 4, "quantile"] <- "4th quartile of gene expression full dataset"
normalized.file.quartiles[normalized.file.quartiles$quantile == 1, "quantile"] <- "1st quartile of gene expression full dataset"

normalized.file <- rbind(normalized.file, normalized.file.quartiles)

# To conserve to same order of groups in the plot : 
ordered.group <- c("random data", "4th quartile of gene expression random data", "1st quartile of gene expression random data", "full dataset", "4th quartile of gene expression full dataset", "1st quartile of gene expression full dataset", "known cycling genes")
normalized.file$quantile <- factor(normalized.file$quantile, levels=ordered.group)

normalized.file.1 <- normalized.file


########################################
# limited.time.points.2 ####################
########################################
file.dir <- paste(paste(main.dir, "DATA", species, tissue, limited.time.points.2, sep = "/"), "/", sep="")
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
### We separate genes en 4 quartiles
raw.species.data <- raw.species.data %>% dplyr::mutate(quantile = ntile(mediane.expr, 4))

normalized.file.quartiles <- merge(normalized.file[, -4], raw.species.data, by = "ID")
# Keep only 1st and 4th quartiles :
normalized.file.quartiles <- subset(normalized.file.quartiles, quantile %in% c(1,4))
normalized.file.quartiles <- normalized.file.quartiles[, c("ID", "pvalue", "algorithm", "quantile")]
normalized.file.quartiles[normalized.file.quartiles$quantile == 4, "quantile"] <- "4th quartile of gene expression full dataset"
normalized.file.quartiles[normalized.file.quartiles$quantile == 1, "quantile"] <- "1st quartile of gene expression full dataset"

normalized.file <- rbind(normalized.file, normalized.file.quartiles)

# To conserve to same order of groups in the plot : 
ordered.group <- c("random data", "4th quartile of gene expression random data", "1st quartile of gene expression random data", "full dataset", "4th quartile of gene expression full dataset", "1st quartile of gene expression full dataset", "known cycling genes")
normalized.file$quantile <- factor(normalized.file$quantile, levels=ordered.group)

normalized.file.2 <- normalized.file

########################################
# limited.time.points.3 ####################
########################################
file.dir <- paste(paste(main.dir, "DATA", species, tissue, limited.time.points.3, sep = "/"), "/", sep="")
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
### We separate genes en 4 quartiles
raw.species.data <- raw.species.data %>% dplyr::mutate(quantile = ntile(mediane.expr, 4))

normalized.file.quartiles <- merge(normalized.file[, -4], raw.species.data, by = "ID")
# Keep only 1st and 4th quartiles :
normalized.file.quartiles <- subset(normalized.file.quartiles, quantile %in% c(1,4))
normalized.file.quartiles <- normalized.file.quartiles[, c("ID", "pvalue", "algorithm", "quantile")]
normalized.file.quartiles[normalized.file.quartiles$quantile == 4, "quantile"] <- "4th quartile of gene expression full dataset"
normalized.file.quartiles[normalized.file.quartiles$quantile == 1, "quantile"] <- "1st quartile of gene expression full dataset"

normalized.file <- rbind(normalized.file, normalized.file.quartiles)

# To conserve to same order of groups in the plot : 
ordered.group <- c("random data", "4th quartile of gene expression random data", "1st quartile of gene expression random data", "full dataset", "4th quartile of gene expression full dataset", "1st quartile of gene expression full dataset", "known cycling genes")
normalized.file$quantile <- factor(normalized.file$quantile, levels=ordered.group)

normalized.file.3 <- normalized.file



colnames(normalized.file.0)[2] <- "pvalue.0"
colnames(normalized.file.1)[2] <- "pvalue.1"
colnames(normalized.file.2)[2] <- "pvalue.2"
colnames(normalized.file.3)[2] <- "pvalue.3"

normalized.file <- normalized.file.0
normalized.file$pvalue.1 <- normalized.file.1$pvalue.1
normalized.file$pvalue.2 <- normalized.file.2$pvalue.2
normalized.file$pvalue.3 <- normalized.file.3$pvalue.3

normalized.file <- subset(normalized.file, quantile == "full dataset")
normalized.file <- subset(normalized.file, algorithm %in% c("ARS", "GeneCycle", "empJTK"))
normalized.file$algorithm <- factor(normalized.file$algorithm, levels=c("ARS", "GeneCycle", "empJTK"))


# Plot
finalPlot <- ggplot(normalized.file, aes(x = pvalue.0, y = pvalue.3)) +
  geom_point(aes(x = -log10(pvalue.0), y = -log10(pvalue.3)), size=0.001, color="grey38") +
  facet_wrap(~ algorithm, ncol = 3) +
  theme_ridges(font_size = 12, grid = TRUE, line_size = 0.1) +
  labs(x = "-log10(default p-value)\n original dataset", y = "-log10(default p-value)\n downsampled dataset") +
  theme(axis.text.x = element_text(size=9, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size=9, hjust = 0.5, vjust = 0.5),
        axis.title.x = element_text(size=10, face="italic", hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(size=10, face="italic", hjust = 0.5, vjust = 0.5),
        strip.text.x = element_text( size=11, face="bold", hjust = 0.5, vjust = 0.5),
        plot.title = element_blank()) 
#
###  Export as pdf file
out.pdf.name <- paste(main.dir, "/RESULTS/", "pvalues_comparison_after_downsampling_for3algo.pdf", sep="")
pdf(file = out.pdf.name, onefile=FALSE, width = 5.1, height = 2.1) # or other device
print(finalPlot)
dev.off()


#### Correlation tests
test <- subset(normalized.file, algorithm == "empJTK")
cor.test(asinh(test$pvalue.0), asinh(test$pvalue.3), na.action="na.exclude", method="pearson") 




