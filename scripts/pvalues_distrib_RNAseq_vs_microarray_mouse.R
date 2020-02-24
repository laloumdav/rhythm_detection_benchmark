########################################
### Density distribution of pvalues ####
############## MOUSE ###################
########## in microarray ###############
############### VS #####################
############# RNAseq ###################
########################################
library(ggplot2)
library(ggridges)
library(cowplot)


args = commandArgs(trailingOnly=TRUE)
tissue <- args[1]

species.microarray <- "mouse_microarray"
species.RNAseq <- "mouse_RNAseq"
species.name <- "mouse"

main.dir <- "~/Documents/rhythm_detection_benchmark"



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
  ########  1) Original microarray dataset     ##################
  ###############################################################
  ###############################################################
  file.dir <- paste(paste(main.dir, "DATA", species.microarray, tissue, sep = "/"), "/", sep="")
  file <- paste(file.dir, "normalized_", p.val, "/normalized_pvalue_per_gene.txt", sep = "")
  normalized.file <- read.table(file, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  # If initial normalization per gene was needed : we keep only the one selected by pvalue.per.gene.normalization.method argument
  if (length(colnames(normalized.file))>3){
    normalized.file <- normalized.file[,c("ID", per.gene.normalization.method, "algorithm")]
  }
  colnames(normalized.file) <- c("ID", "pvalue", "algorithm")
  normalized.file$dataset <- rep("Microarray (original dataset: each 2h over 48h)", nrow(normalized.file))
  

  ###############################################################
  ###############################################################
  ########  2) Restricted time-points microarray data     #######
  ###############################################################
  ###############################################################
  file.dir <-  paste(file.dir, "exploration_within_restricted_timepoints/", sep="")
  file <- paste(file.dir, "normalized_", p.val, "/normalized_pvalue_per_gene.txt", sep = "")
  normalized.file.bis <- read.table(file, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  # If initial normalization per gene was needed : we keep only the one selected by pvalue.per.gene.normalization.method argument
  if (length(colnames(normalized.file.bis))>3){
    normalized.file.bis <- normalized.file.bis[,c("ID", per.gene.normalization.method, "algorithm")]
  }
  colnames(normalized.file.bis) <- c("ID", "pvalue", "algorithm")
  normalized.file.bis$dataset <- rep("Microarray (restricted timepoints: each 6h over 48h)", nrow(normalized.file.bis))
  
  
  ###########################
  # Final Table
  normalized.file <- rbind(normalized.file, normalized.file.bis)
  
  
  
  
  ###############################################################
  ###############################################################
  ########  3) Original RNAseq dataset     ######################
  ###############################################################
  ###############################################################
  file.dir <- paste(paste(main.dir, "DATA", species.RNAseq, tissue, sep = "/"), "/", sep="")
  
  file <- paste(file.dir, "normalized_", p.val, "/normalized_pvalue_per_gene.txt", sep = "")
  normalized.file.bis <- read.table(file, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  # If initial normalization per gene was needed : we keep only the one selected by pvalue.per.gene.normalization.method argument
  if (length(colnames(normalized.file.bis))>3){
    normalized.file.bis <- normalized.file.bis[,c("ID", per.gene.normalization.method, "algorithm")]
  }
  colnames(normalized.file.bis) <- c("ID", "pvalue", "algorithm")
  normalized.file.bis$dataset <- rep("RNAseq (original dataset: each 6h over 48h)", nrow(normalized.file.bis))
  
  ###########################
  # Final Table
  normalized.file <- rbind(normalized.file, normalized.file.bis)
  
  ###############################################################
  ###############################################################
  ########  4) RNAseq data with low expressed genes removed #####
  ###############################################################
  ###############################################################
  #file.dir <- paste(paste(main.dir, "DATA", "mouse_RNAseq", tissue, sep = "/"), "/", sep="")
  #file.dir <-  paste(file.dir, "exploration_through_dataset_with_low_expressed_genes_removed/", sep="")
  #file <- paste(file.dir, "normalized_", p.val, "/normalized_pvalue_per_gene.txt", sep = "")
  #normalized.file.bis <- read.table(file, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  # If initial normalization per gene was needed : we keep only the one selected by pvalue.per.gene.normalization.method argument
  #if (length(colnames(normalized.file.bis))>3){
  #  normalized.file.bis <- normalized.file.bis[,c("ID", per.gene.normalization.method, "algorithm")]
  #}
  #colnames(normalized.file.bis) <- c("ID", "pvalue", "algorithm")
  #normalized.file.bis$dataset <- rep("RNAseq_with_low_expressed_genes_removed", nrow(normalized.file.bis))
  
  ###########################
  # Final Table
  #normalized.file <- rbind(normalized.file, normalized.file.bis)
  
  
  # To conserve to same order of groups in the plot : 
  ordered.group <- c("Microarray (original dataset: each 2h over 48h)", "Microarray (restricted timepoints: each 6h over 48h)", "RNAseq (original dataset: each 6h over 48h)")
  normalized.file$dataset <- factor(normalized.file$dataset, levels=ordered.group)
  # Re-order the order of algorithms
  ordered.group <- c("ARS", "GeneCycle", "JTK", "LS", "RAIN", "empJTK", "meta2d")
  normalized.file$algorithm <- factor(normalized.file$algorithm, levels=ordered.group)
  
  
  pval.text1 <- gsub(".pvalue", "\n p-value", p.val)
  pval.text2 <- gsub(".pvalue", " p-value", p.val)
  ##########################
  ######### PLOT ###########
  ##########################
  densityDistribPlot <- ggplot(normalized.file, aes(x = pvalue, fill = fct_rev(factor(dataset)), color = NA)) +
    geom_density(aes(x = pvalue, fill = factor(dataset)), alpha=0.7, size=0.2, color="black") +
    scale_fill_manual(values = c(`Microarray (original dataset: each 2h over 48h)`="wheat3", `Microarray (restricted timepoints: each 6h over 48h)`="sandybrown", `RNAseq (original dataset: each 6h over 48h)`="steelblue3")) +
    facet_wrap(~ algorithm, scales = "free", ncol = 2) +
    scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), limits = c(0, 1)) +
    labs(title = "Density distribution",
         x = pval.text2, y = "density") +
    theme_ridges(grid = TRUE, line_size = 0.1) +
    theme(axis.text.x = element_text(size=15, face="bold.italic", hjust = 0.5, vjust = 0.5),
          axis.title.x = element_text(size=16, face="bold.italic", hjust = 0.5, vjust = 0.5),
          axis.title.y = element_text(size=16, face="bold", hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size=15, hjust = 0.5, vjust = 0.5),
          strip.text = element_text(size=15, face="bold.italic", hjust = 0.5, vjust = 0.5),
          strip.background = element_rect(fill="transparent", color=NA),
          plot.title = element_text(size=13, face="bold", hjust = 0.5),
          legend.position = c(0.5, 0.14),
          legend.direction = 'vertical',
          legend.title = element_blank(),
          legend.text = element_text(size=13, face="bold.italic"),
          legend.key = element_rect(color = "transparent", fill = "transparent"),
          legend.spacing.x = unit(0.3, "cm"),
          legend.spacing.y = unit(0.3, "cm"),
          axis.ticks.y = element_blank())
  
  
  
  
  
  
  ###  Export as pdf file
  logo <- paste("~/Documents/workspace/DATA/images/", species.name, "_image.png", sep="")
  out.pdf.name <- paste(main.dir, "/RESULTS/microarray_vs_RNAseq/", p.val, "_distrib_", 
                        species.name, "_", tissue, ".pdf", sep="")
  tissue.text <- gsub("_", " ", tissue)
  final.draw.plot <- ggdraw() +
    draw_plot(densityDistribPlot) +
    draw_image(logo, #scale = 0.2, 
               x = 0.8, y = 0.008,
               width = 0.15, height = 0.15) +
    draw_label(tissue.text, colour = "black", fontface = "bold", size = 10, 
               x = 0.88, y = 0.02)
  
  pdf(file = out.pdf.name, onefile=FALSE, width = 9, height = 10) # or other device
  print(final.draw.plot)
  dev.off()
  
}
