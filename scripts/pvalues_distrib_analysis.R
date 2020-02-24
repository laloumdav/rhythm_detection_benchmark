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


pvalues <- c("raw.pvalue", "default.pvalue", "BH.Q")
for (a in 1:length(pvalues)) {
  
  p.val <- pvalues[a]
  
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
  ### We separate genes en 4 quartiles
  raw.species.data <- raw.species.data %>% dplyr::mutate(quantile = ntile(mediane.expr, 4))
  
  normalized.file.quartiles <- merge(normalized.file[, -4], raw.species.data, by = "ID")
  # Keep only 1st and 4th quartiles :
  normalized.file.quartiles <- subset(normalized.file.quartiles, quantile %in% c(1,4))
  normalized.file.quartiles <- normalized.file.quartiles[, c("ID", "pvalue", "algorithm", "quantile")]
  normalized.file.quartiles[normalized.file.quartiles$quantile == 4, "quantile"] <- "4th quartile of gene expression full dataset"
  normalized.file.quartiles[normalized.file.quartiles$quantile == 1, "quantile"] <- "1st quartile of gene expression full dataset"
  
  
  normalized.file <- rbind(normalized.file, normalized.file.quartiles)
  
  

  
  ###############################################################
  ###############################################################
  ##########  2) On random data     #############################
  ###############################################################
  ###############################################################
  
  file.dir <- paste(paste(main.dir.ramdom.data, "DATA", species, tissue, sep = "/"), "/", sep="")
  
  file <- paste(file.dir, "normalized_", p.val, "/normalized_pvalue_per_gene.txt", sep = "")
  normalized.file.bis <- read.table(file, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  # If initial normalization per gene was needed : we keep only the one selected by pvalue.per.gene.normalization.method argument
  if (length(colnames(normalized.file.bis))>3){
    normalized.file.bis <- normalized.file.bis[,c("ID", per.gene.normalization.method, "algorithm")]
  }
  colnames(normalized.file.bis) <- c("ID", "pvalue", "algorithm")
  normalized.file.bis$quantile <- rep("random data", nrow(normalized.file.bis))
  
  
  ###########################
  # Final Table
  normalized.file <- rbind(normalized.file, normalized.file.bis)
  
  
  
  
  
  
  ####################
  # Through quartiles of RANDOM DATA : 
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
  
  raw.species.data[raw.species.data$quantile == 4, "quantile"] <- "4th quartile of gene expression random data"
  raw.species.data[raw.species.data$quantile == 1, "quantile"] <- "1st quartile of gene expression random data"
  
  normalized.file.quartiles <- merge(normalized.file.bis[, -4], raw.species.data[, -2], by = "ID")
  # Keep only 1st and 4th quartiles :
  normalized.file.quartiles <- subset(normalized.file.quartiles, quantile %in% c("4th quartile of gene expression random data", "1st quartile of gene expression random data"))
  normalized.file.quartiles <- normalized.file.quartiles[, c("ID", "pvalue", "algorithm", "quantile")]
  
  
  
  normalized.file <- rbind(normalized.file, normalized.file.quartiles)
  
  
  
  
  
  
  
  ###############################################################
  ###############################################################
  ######  3) Through known circad genes   #######################
  ######    (if the genes list exist)     #######################
  ###############################################################
  ###############################################################
  file.dir <- paste(paste(main.dir, "DATA", species, tissue, sep = "/"), "/", sep="")
  file.dir <-  paste(file.dir, "exploration_through_known_circad_genes/", sep="")
  file <- paste(file.dir, "normalized_", p.val, "/normalized_pvalue_per_gene.txt", sep = "")
  
  if (file.exists(file)) {
    normalized.file.bis <- read.table(file, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    # If initial normalization per gene was needed : we keep only the one selected by pvalue.per.gene.normalization.method argument
    if (length(colnames(normalized.file.bis))>3){
      normalized.file.bis <- normalized.file.bis[,c("ID", per.gene.normalization.method, "algorithm")]
    }
    colnames(normalized.file.bis) <- c("ID", "pvalue", "algorithm")
    normalized.file.bis$quantile <- rep("known cycling genes", nrow(normalized.file.bis))
    
    
    ###########################
    # Final Table
    normalized.file <- rbind(normalized.file, normalized.file.bis) 
  }

  # To conserve to same order of groups in the plot : 
  ordered.group <- c("random data", "4th quartile of gene expression random data", "1st quartile of gene expression random data", "full dataset", "4th quartile of gene expression full dataset", "1st quartile of gene expression full dataset", "known cycling genes")
  normalized.file$quantile <- factor(normalized.file$quantile, levels=ordered.group)
  
  # Re-order the order of algorithms
  ordered.group <- c("ARS", "GeneCycle", "JTK", "LS", "RAIN", "empJTK", "meta2d")
  normalized.file$algorithm <- factor(normalized.file$algorithm, levels=ordered.group)
  
  

  pval.text <- gsub(".pvalue", " p-value", p.val)
  # Plot
  pValuesDistributionPlot <-  ggplot(normalized.file, aes(x = pvalue, y = fct_rev(quantile), fill = quantile, color=NA)) +
    geom_density_ridges_gradient(scale=1.8) +
    facet_wrap(~ algorithm, scales = "free", ncol = 2) +
    scale_fill_manual(values = c(`random data`="wheat3", `4th quartile of gene expression full dataset`="aquamarine4", `1st quartile of gene expression full dataset`="darkseagreen", 
                                 `4th quartile of gene expression random data`="sienna3", `1st quartile of gene expression random data`="sienna1", `full dataset`="grey38", `known cycling genes`="indianred")) +
    scale_color_manual(values = c(`random data`="wheat3", `4th quartile of gene expression full dataset`="aquamarine4", `1st quartile of gene expression full dataset`="darkseagreen", 
                                `4th quartile of gene expression random data`="sienna3", `1st quartile of gene expression random data`="sienna1", `full dataset`="grey38", `known cycling genes`="indianred")) +
    theme_ridges(font_size = 12, grid = TRUE, line_size = 0.1) +
    labs(title = "Density distribution", x = pval.text, y = "density") +
    scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), limits = c(0, 1)) +
    theme(axis.text.x = element_text(size=9, face="bold.italic", hjust = 0.5, vjust = 0.5),
          axis.title.x = element_text(size=10, face="bold.italic", hjust = 0.5, vjust = 0.5),
          axis.title.y = element_text(size=10, face="bold", hjust = 0.5, vjust = 0.5),
          strip.text.x = element_text( size=11, face="bold", hjust = 0.5, vjust = 0.5),
          strip.background = element_rect(fill = "transparent", color=NA),
          plot.title = element_blank(),#element_text(size=13, face="bold", hjust = 0.5),
          legend.position = c(0.5, 0.15),
          legend.direction = 'vertical',
          legend.title = element_blank(),
          legend.text = element_text( size=8, face="bold", hjust = 0),
          legend.spacing.x = unit(0.2, "cm"),
          legend.spacing.y = unit(0.4, "cm"),
          legend.key.height = unit(0.4, "cm"),
          #axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) 



  
  
  ###  Export as pdf file
  logo <- paste(main.dir, "/DATA/images/", species.name, "_image.png", sep="")
  out.pdf.name <- paste(main.dir, "/RESULTS/pval_distrib_images/", p.val, "_distrib_", 
                        species, "_", tissue, ".pdf", sep="")
  tissue.text <- gsub("_", " ", tissue)
  final.draw.plot <- ggdraw() +
    draw_plot(pValuesDistributionPlot) +
    draw_image(logo, #scale = 0.2, 
               x = 0.8, y = 0.008,
               width = 0.15, height = 0.15) +
    draw_label(tissue.text, colour = "black", fontface = "bold", size = 10, 
               x = 0.88, y = 0.02)
  
  pdf(file = out.pdf.name, onefile=FALSE, width = 6, height = 8) # or other device
  print(final.draw.plot)
  dev.off()
  
  # Also save it as a R object:
  out.rds.name <- paste(main.dir, "/RESULTS/pval_distrib_images/", p.val, "_distrib_", 
                        species, "_", tissue, ".rds", sep="")
  saveRDS(final.draw.plot, file = out.rds.name)
  
  
  
}



