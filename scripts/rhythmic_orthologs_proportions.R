##############################
##############################
###### proportion PLOTS ######
##############################
##############################
library(reshape2)
library(ggplot2)
require(RColorBrewer)
library(tidyverse)
library(wesanderson)
library(ggridges)
library(cowplot)
library(magick)
library(ggforce)


args = commandArgs(trailingOnly=TRUE)
species1 <- args[1]
tissue1 <- args[2]
species2 <- args[3]
tissue2 <- args[4]
species2.algorithm <- args[5]
species2.threshold = args[6]
species2.pvalue <- args[7]
point.species1.threshold = args[8]
only.within.conserved.genes <- as.logical(args[9])
only.ortho.one.to.one <- as.logical(args[10])

point.species1.FDR.threshold <- FALSE

main.dir <- "~/Documents/rhythm_detection_benchmark"


if (grepl("mouse", species1)){ species1.name <- "mouse" } else { species1.name <- species1 }
if (grepl("mouse", species2)){ species2.name <- "mouse" } else { species2.name <- species2 }


p.val <- species2.pvalue

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





########################################################################
# Create a "Naive" method to be compared with others 
# Which will order genes (or just orthologs) according to their mean expression 
# From high expressed to low expressed gene
# And for each gene i, will calculate the proportion of rhythmic orthologs within genes set (from high expressed to this gene i)
########################################################################
# Retreive original data
raw.species1.file.name <- paste(main.dir, "/DATA/", species1, "/", tissue1, "/", tissue1, ".txt", sep = "")
raw.species1.data <- read.table(raw.species1.file.name, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
# Keep only genes available in master object
raw.species1.data <- raw.species1.data[which(raw.species1.data$ID %in% master[, paste(species1.name, "ID", sep="_")]), ]
# calculate the mean of expression for each gene
raw.species1.data$median.expr <- apply(raw.species1.data[, -1], 1, min)
raw.species1.data <- raw.species1.data[, c("ID", "median.expr")]
# if several transcripts per gene, we keep the mean of mediane gene expr to get 1 value per gene
raw.species1.data <- aggregate(raw.species1.data["median.expr"], raw.species1.data[1], max)
colnames(raw.species1.data)[1] <- paste(species1.name, "ID", sep="_")
naive.method <- merge(raw.species1.data, master, by=paste(species1.name, "ID", sep="_"))
if (only.within.conserved.genes) {
  naive.method <- subset(naive.method, conserved_gene == TRUE)
}
if (only.ortho.one.to.one) {
  naive.method <- subset(naive.method, orthotype == "1:1")
}
# Order genes from high expressed to low expressed genes
naive.method <- naive.method[order(naive.method[, "median.expr"], decreasing = TRUE), ]
naive.method$genes_ordered <- 1:nrow(naive.method)


# orthologs list 
orthologs.list <- naive.method[, paste(species1.name, "ID", sep="_")]
# rhythmic orthologs list (orthologs rhythmic in species_2)
rhythmic.orthologs.list <- naive.method[naive.method[[paste(species2.name, species2.algorithm, sep="_")]] <= species2.threshold, paste(species1.name, "ID", sep="_")]

naive.method <- naive.method[, c(paste(species1.name, "ID", sep="_"), "genes_ordered")]
for (i in 1:nrow(naive.method)) {
  genes.list <- naive.method[1:i, 1]
  rhythmic.orthologs.nb <- length(intersect(genes.list, rhythmic.orthologs.list))
  naive.method$proportion[i] <- rhythmic.orthologs.nb/i
}
naive.method$algorithm <- "Naive"
naive.method[[species2.pvalue]] <- NA



###
### Tot number of orthologs 
###
sub.master <- master
if (only.within.conserved.genes) {
  sub.master <- subset(sub.master, conserved_gene == TRUE)
}
if (only.ortho.one.to.one) {
  sub.master <- subset(sub.master, orthotype == "1:1")
}
tot.orthologs.nb = length(orthologs.list)


###
### For each algorithm, we order genes according to their p-value 
### Then, for each ordered gene x, we count the number of genes with a p-value < p-value of the gene x  = A
### And we count the number of these genes which are detected rhythmic in species2  = B
### And calculate the ratio B/A
### 
ordered.genes <- 1:nrow(sub.master)
ratio.table <- data.frame(genes_ordered = ordered.genes)
pval.table <- data.frame(genes_ordered = ordered.genes)
for (j in 1:length(algorithms.list)) {
  sub.master <- sub.master[order(sub.master[, algorithms.list[j]]), ]
  
  for (x in 1:nrow(sub.master)) {
    pvals.species2 <- sub.master[1:x, paste(species2.name, species2.algorithm, sep="_")]
    nb.rhythmic.species2 <- length(which(pvals.species2 <= species2.threshold))
    sub.master[x, "ratio"] <- nb.rhythmic.species2/x
  }
  
  ratio.table[, algorithms.list[j]] <- sub.master[, "ratio"]
  pval.table[, algorithms.list[j]] <- sub.master[, algorithms.list[j]]
}

proportions.table.final <- cbind(stack(ratio.table, select=algorithms.list),
                                 stack(pval.table, select=algorithms.list))[, c(1,3,4)]
colnames(proportions.table.final) <- c("proportion", p.val, "algorithm") 
proportions.table.final$genes_ordered <- rep(ordered.genes, length(algorithms.list))

colnames(naive.method)[5] <- p.val
naive.method <- naive.method[, colnames(proportions.table.final)]
proportions.table.final <- rbind(proportions.table.final, naive.method)



# To conserve to same order of groups in the plot : 
ordered.group <- c("ARS", "GeneCycle", "empJTK", "JTK", "LS", "meta2d", "RAIN", "Naive")
proportions.table.final$algorithm <- factor(proportions.table.final$algorithm, levels=ordered.group)


##################
###### PLOT ######
################## 
colors.list <- c(wes_palette("Darjeeling1"), wes_palette("GrandBudapest2"))

proportionPlot <- ggplot(proportions.table.final, aes(x=conserved_genes_nb_for_ideal.cutoff, y=proportion_for_ideal.cutoff)) +
  #geom_line(aes(x=genes_ordered, y=proportion, col=algorithm), size=0.8) +
  geom_smooth(aes(x=genes_ordered, y=proportion, col=algorithm), method = "loess", span = 0.02, se = FALSE, size=0.7) +
  #geom_point(data = best.cutoff, aes(col=algorithm, fill=algorithm), na.rm = TRUE, 
  #           size = 4, shape=23, color="black") +
  theme_light(base_size = 0.25) +
  ylim(0.2, 0.9) +
  scale_fill_manual(values = c(JTK=colors.list[1], LS=colors.list[6], empJTK=colors.list[7], ARS=colors.list[4], RAIN=colors.list[5], meta2d=colors.list[2], GeneCycle="blue", Naive="black")) +
  scale_color_manual(values = c(JTK=colors.list[1], LS=colors.list[6], empJTK=colors.list[7], ARS=colors.list[4], RAIN=colors.list[5], meta2d=colors.list[2], GeneCycle="blue", Naive="black")) +
  guides(colour = guide_legend(override.aes = list(shape = NA, size = 2))) +
  labs(x = "number of orthologs detected rhythmic", y = paste("proportion which are\nrhythmic in", species2.name, tissue2, sep=" ")) +
  theme(axis.text.x = element_text(size=12, hjust = 0.5, vjust = 0.5, margin = margin(l = 20)),
        axis.title.x = element_text(size=14, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size=12, hjust = 0.5, vjust = 0.5, margin = margin(l = 20)),
        axis.title.y = element_text(size=14, hjust = 0.5, vjust = 0.5),
        legend.position = c(0.87, 0.75), 
        legend.title = element_text(size=13, face="bold.italic", colour="darkorchid4", hjust = 0.5, vjust = 0.5),
        legend.text = element_text(size=12, hjust = 0, vjust = 0.5))

#


###
### Proportion plot with dots corresponding to usual cutoffs for each algorithm
###
usual.cutoff <- subset(proportions.table.final, proportions.table.final[[p.val]] <= point.species1.threshold)

cutoff.table <- data.frame(algorithm=algorithms.list)
for (i in 1:length(algorithms.list)) {
  ortho.nb.cutoff <- max(subset(usual.cutoff, algorithm == algorithms.list[i])$genes_ordered)
  if (!is.infinite(ortho.nb.cutoff)) {
    proportion.cutoff <- subset(usual.cutoff, algorithm == algorithms.list[i] & genes_ordered==ortho.nb.cutoff)$proportion
    cutoff.table[cutoff.table$algorithm==algorithms.list[i], "cutoff"] <- point.species1.threshold
    cutoff.table[cutoff.table$algorithm==algorithms.list[i], "ortho.nb.cutoff"] <- ortho.nb.cutoff
    cutoff.table[cutoff.table$algorithm==algorithms.list[i], "proportion.cutoff"] <- proportion.cutoff
  }
}


### --- FDR step --- ###
#point.species1.FDR.threshold = TRUE # Remove "#" to get FDR thresholds in the plot
#FDR.cutoff.1 = 0.005 # Remove "#" to get FDR thresholds in the plot
#FDR.cutoff.2 = 0.05 # Remove "#" to get FDR thresholds in the plot

if (point.species1.FDR.threshold == TRUE) {

  proportions.table.final$FDR <- p.adjust(proportions.table.final[[p.val]], method = "fdr")
  
  usual.cutoff <- subset(proportions.table.final, FDR <= FDR.cutoff.1)
  cutoff.table.bis <- data.frame(algorithm=algorithms.list)
  for (i in 1:length(algorithms.list)) {
    ortho.nb.cutoff <- max(subset(usual.cutoff, algorithm == algorithms.list[i])$genes_ordered)
    if (!is.infinite(ortho.nb.cutoff)) {
      proportion.cutoff <- subset(usual.cutoff, algorithm == algorithms.list[i] & genes_ordered==ortho.nb.cutoff)$proportion
      cutoff.table.bis[cutoff.table.bis$algorithm==algorithms.list[i], "cutoff"] <- point.species1.threshold
      cutoff.table.bis[cutoff.table.bis$algorithm==algorithms.list[i], "ortho.nb.cutoff"] <- ortho.nb.cutoff
      cutoff.table.bis[cutoff.table.bis$algorithm==algorithms.list[i], "proportion.cutoff"] <- proportion.cutoff
    }
  }
  
  usual.cutoff <- subset(proportions.table.final, FDR <= FDR.cutoff.2)
  cutoff.table.bis.2 <- data.frame(algorithm=algorithms.list)
  for (i in 1:length(algorithms.list)) {
    ortho.nb.cutoff <- max(subset(usual.cutoff, algorithm == algorithms.list[i])$genes_ordered)
    if (!is.infinite(ortho.nb.cutoff)) {
      proportion.cutoff <- subset(usual.cutoff, algorithm == algorithms.list[i] & genes_ordered==ortho.nb.cutoff)$proportion
      cutoff.table.bis.2[cutoff.table.bis.2$algorithm==algorithms.list[i], "cutoff"] <- point.species1.threshold
      cutoff.table.bis.2[cutoff.table.bis.2$algorithm==algorithms.list[i], "ortho.nb.cutoff"] <- ortho.nb.cutoff
      cutoff.table.bis.2[cutoff.table.bis.2$algorithm==algorithms.list[i], "proportion.cutoff"] <- proportion.cutoff
    }
  }
  
}
### --- FDR step --- ###


###
### Proportion plot with dots corresponding to ideal cutoffs for each algorithm
###
proportionPlotWithDots <- ggplot(proportions.table.final, aes(x=ortho.nb.cutoff, y=proportion.cutoff)) +
  #geom_line(aes(x=genes_ordered, y=proportion, col=algorithm), size=0.8) +
  geom_smooth(aes(x=genes_ordered, y=proportion, col=algorithm), method = "loess", span = 0.06, se = FALSE, size=0.7) +
  geom_point(data = cutoff.table, aes(col=algorithm, fill=algorithm), na.rm = TRUE, 
             size = 2.5, shape=23, color="black") +
  theme_light(base_size = 0.25) +
  scale_fill_manual(values = c(JTK=colors.list[1], LS=colors.list[6], empJTK=colors.list[7], ARS=colors.list[4], RAIN=colors.list[5], meta2d=colors.list[2], GeneCycle="blue", Naive="black")) +
  scale_color_manual(values = c(JTK=colors.list[1], LS=colors.list[6], empJTK=colors.list[7], ARS=colors.list[4], RAIN=colors.list[5], meta2d=colors.list[2], GeneCycle="blue", Naive="black")) +
  guides(colour = guide_legend(override.aes = list(shape = NA, size = 2))) +
  labs(x = "x\nnumber of orthologs detected rhythmic", y = paste("proportion which are\nrhythmic in", species2.name, tissue2, sep=" ")) +
  theme(axis.text.x = element_text(size=12, hjust = 0.5, vjust = 0.5, margin = margin(l = 20)),
        axis.title.x = element_text(size=14, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size=12, hjust = 0.5, vjust = 0.5, margin = margin(l = 20)),
        axis.title.y = element_text(size=14, hjust = 0.5, vjust = 0.5),
        legend.position = "none")


if (point.species1.FDR.threshold == TRUE) {
  proportionPlotWithDots <- ggplot(proportions.table.final, aes(x=ortho.nb.cutoff, y=proportion.cutoff)) +
    #geom_line(aes(x=genes_ordered, y=proportion, col=algorithm), size=0.8) +
    geom_smooth(aes(x=genes_ordered, y=proportion, col=algorithm), method = "loess", span = 0.06, se = FALSE, size=0.7) +
    geom_point(data = cutoff.table, aes(col=algorithm, fill=algorithm), na.rm = TRUE, 
               size = 3, shape=21, color="black") +
    geom_point(data = cutoff.table.bis, aes(col=algorithm, fill=algorithm), na.rm = TRUE, 
               size = 2.5, shape=22, color="black") +
    geom_point(data = cutoff.table.bis.2, aes(col=algorithm, fill=algorithm), na.rm = TRUE, 
               size = 3, shape=23, color="black") +
    #ylim(0.1, 0.9) +
    theme_light(base_size = 0.25) +
    scale_fill_manual(values = c(JTK=colors.list[1], LS=colors.list[6], empJTK=colors.list[7], ARS=colors.list[4], RAIN=colors.list[5], meta2d=colors.list[2], GeneCycle="blue", Naive="black")) +
    scale_color_manual(values = c(JTK=colors.list[1], LS=colors.list[6], empJTK=colors.list[7], ARS=colors.list[4], RAIN=colors.list[5], meta2d=colors.list[2], GeneCycle="blue", Naive="black")) +
    guides(colour = guide_legend(override.aes = list(shape = NA, size = 2))) +
    labs(x = "x\nnumber of orthologs detected rhythmic", y = paste("proportion which are\nrhythmic in", species2.name, tissue2, sep=" ")) +
    theme(axis.text.x = element_text(size=12, hjust = 0.5, vjust = 0.5, margin = margin(l = 20)),
          axis.title.x = element_text(size=14, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size=12, hjust = 0.5, vjust = 0.5, margin = margin(l = 20)),
          axis.title.y = element_text(size=14, hjust = 0.5, vjust = 0.5),
          legend.position = "none")
}

#


###  Export as pdf file
pval.text1 <- "x = f(cutoff)"
pval.text2 <- gsub(".pvalue", "\np-value", species2.pvalue)
pval.text2 <- paste(pval.text2, "<", species2.threshold, sep="")
ratio.illust <- paste(main.dir, "/DATA/images/ratio_illustration.png", sep="")
logo1 <- paste(main.dir, "/DATA/images/", species1.name, "_image.png", sep="")
logo2 <- paste(main.dir, "/DATA/images/", species2.name, "_image.png", sep="")
legend <- paste(main.dir, "/DATA/images/ratio_legend_bis.png", sep="")
usual.cutoff.text <- paste("cutoff", point.species1.threshold , sep=" = ")

out.pdf.name <- paste(main.dir, "/RESULTS/proportion_rhythmic_orthologs_images/", species2.algorithm, "_", species2.pvalue, "_RATIO_", 
                      species1, "_", species2, "_", tissue1, ".pdf", sep="")
final.draw.plot <- ggdraw() +
  draw_plot(proportionPlotWithDots, 
            x = 0.05, y = 0.05, width = 0.7, height = 0.6) +
  
  draw_image(ratio.illust, scale = 1, 
             x = -0.13, y = 0.7, width = 0.77, height = 0.2) +
  draw_image(logo1, scale = 1, 
             x = 0.138, y = 0.77, width = 0.1, height = 0.1) +
  draw_image(logo2, scale = 1, 
             x = 0.25, y = 0.77, width = 0.1, height = 0.1) +
  draw_label("method", colour = "darkorchid4", fontface = "bold.italic", size = 11, 
             x = 0.12, y = 0.86) +
  draw_label(pval.text1, colour = "grey39", fontface = "italic", size = 11, 
             x = 0.125, y = 0.76) +
  draw_label(species2.algorithm, colour = "black", fontface = "bold", size = 11, 
             x = 0.39, y = 0.86) +
  draw_label(pval.text2, colour = "black", fontface = "bold", size = 10, 
             x = 0.38, y = 0.76) +
  draw_label(tissue1, colour = "black", fontface = "bold", size = 11, 
             x = 0.21, y = 0.72) +
  draw_label(tissue2, colour = "black", fontface = "bold", size = 11, 
             x = 0.3, y = 0.72) +
  
  draw_image(legend, scale = 1.5, 
             x = 0.5, y = 0.7, width = 0.5, height = 0.2) +
  draw_label(usual.cutoff.text, colour = "grey39", fontface = "italic", size = 12, 
             x = 0.755, y = 0.675) 
  
pdf(file = out.pdf.name, onefile=FALSE, width = 11, height = 10) # or other device
print(final.draw.plot)
dev.off()

# Also save it as a R object for Rmarkdown:
out.rds.name <- paste(main.dir, "/RESULTS/proportion_rhythmic_orthologs_images/", species2.algorithm, "_", species2.pvalue, "_RATIO_", 
                      species1, "_", species2, "_", tissue1, ".rds", sep="")
saveRDS(final.draw.plot, file = out.rds.name)

