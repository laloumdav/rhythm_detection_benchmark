######################################################################
######################################################################
# log(meadian(gene expression)) VS p-values given by algorithms ######
########## Comparing tissues #######################################
######################################################################
######################################################################
library(tidyverse)
library(ggplot2)
require(RColorBrewer)
library(reshape2)
library(magick)
library(cowplot)

main.dir <- "~/Documents/rhythm_detection_benchmark"



######
species <- "mouse_microarray"
tissue.list <- c("liver", "lung", "kidney", "muscle", "aorta")
p.val <- "default.pvalue"
######
if (grepl("mouse", species)){ species.name <- "mouse" } else { species.name <- species }



full.dataframe <- data.frame()
for (t in 1:length(tissue.list)) {
  tissue <- tissue.list[t]
  
  file.dir <- paste(paste(main.dir, "DATA", species, tissue, sep = "/"), "/", sep="")
  
  algorithms.list <- list.files(file.dir)
  algorithms.list <- algorithms.list[ grepl("GeneCycle|empJTK|RAIN|JTK|LS|ARS|meta2d", algorithms.list) ]
  algorithms.list.name <- gsub(".txt", "", algorithms.list)
  algorithms.list.files <- paste(file.dir, algorithms.list, sep = "")
  
  file.name <- paste(file.dir, tissue, ".txt", sep = "")
  raw.data <- read.table(file.name, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  raw.data$median <- apply(raw.data[, grep("CT|ZT", colnames(raw.data))], 1, FUN = function(x) round(median(x), digits = 1))
  
  algorithm <- NULL
  for (i in 1:length(algorithms.list)){
    algorithm[[i]] <- read.table(algorithms.list.files[i], head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    algorithm[[i]]$algorithm <- rep(algorithms.list.name[i], nrow(algorithm[[i]]))
    algorithm[[i]]$expr_value <- raw.data$median
  }
  # We make a unique dataframe with these data:
  full.dataframe.tmp <- data.frame()
  for (i in 1:length(algorithms.list)){
    full.dataframe.tmp <- rbind(full.dataframe.tmp, algorithm[[i]])
  }
  
  full.dataframe.tmp$tissue <- rep(tissue, nrow(full.dataframe.tmp))
  full.dataframe <- rbind(full.dataframe, full.dataframe.tmp)

}


# To conserve to same order of groups in the plot : 
ordered.group <- tissue.list
full.dataframe$tissue <- factor(full.dataframe$tissue, levels=ordered.group)


algoTissueSubset.full.dataframe.ARS <- list()
algoTissueSubset.full.dataframe.empJTK <- list()
algoTissueSubset.full.dataframe.GeneCycle <- list()
algoTissueSubset.full.dataframe.JTK <- list()
algoTissueSubset.full.dataframe.LS <- list()
algoTissueSubset.full.dataframe.meta2d <- list()
algoTissueSubset.full.dataframe.RAIN <- list()

tissues.list <- c("aorta", "kidney", "liver", "lung", "muscle")
for (i in 1:length(tissues.list)) {
  tissueSubset.full.dataframe <- subset(full.dataframe, tissue == tissues.list[i])
  algoTissueSubset.full.dataframe.ARS[[i]] <- subset(tissueSubset.full.dataframe, algorithm == "ARS")
  algoTissueSubset.full.dataframe.empJTK[[i]] <- subset(tissueSubset.full.dataframe, algorithm == "empJTK")
  algoTissueSubset.full.dataframe.GeneCycle[[i]] <- subset(tissueSubset.full.dataframe, algorithm == "GeneCycle")
  algoTissueSubset.full.dataframe.JTK[[i]] <- subset(tissueSubset.full.dataframe, algorithm == "JTK")
  algoTissueSubset.full.dataframe.LS[[i]] <- subset(tissueSubset.full.dataframe, algorithm == "LS")
  algoTissueSubset.full.dataframe.meta2d[[i]] <- subset(tissueSubset.full.dataframe, algorithm == "meta2d")
  algoTissueSubset.full.dataframe.RAIN[[i]] <- subset(tissueSubset.full.dataframe, algorithm == "RAIN")
}



########
# PLOT #
########
out.pdf.name <- paste(main.dir, "/RESULTS/geneExpr_vs_pval/", "geneExpr_vs_pval_ARS_", 1, ".pdf", sep="")
pdf(file = out.pdf.name, onefile=FALSE, width = 3, height = 3) # or other device
par(mfcol = c(1, 1), mar=c(0.1, 2.2, 0.1, 0.1))
smoothScatter(log2(algoTissueSubset.full.dataframe.ARS[[1]]$expr_value), -log10(algoTissueSubset.full.dataframe.ARS[[1]]$default.pvalue), xlab = NA, ylab=NA, cex = 2)
dev.off()

out.pdf.name <- paste(main.dir, "/RESULTS/geneExpr_vs_pval/", "geneExpr_vs_pval_empJTK_", 1, ".pdf", sep="")
pdf(file = out.pdf.name, onefile=FALSE, width = 3, height = 3) # or other device
par(mfcol = c(1, 1), mar=c(0.1, 2.2, 0.1, 0.1))
smoothScatter(log2(algoTissueSubset.full.dataframe.empJTK[[1]]$expr_value), -log10(algoTissueSubset.full.dataframe.empJTK[[1]]$default.pvalue), xlab = NA, ylab=NA,cex = 2)
dev.off()

out.pdf.name <- paste(main.dir, "/RESULTS/geneExpr_vs_pval/", "geneExpr_vs_pval_GeneCycle_", 1, ".pdf", sep="")
pdf(file = out.pdf.name, onefile=FALSE, width = 3, height = 3) # or other device
par(mfcol = c(1, 1), mar=c(0.1, 2.2, 0.1, 0.1))
smoothScatter(log2(algoTissueSubset.full.dataframe.GeneCycle[[1]]$expr_value), -log10(algoTissueSubset.full.dataframe.GeneCycle[[1]]$default.pvalue), xlab = NA, ylab=NA,cex = 2)
dev.off()

out.pdf.name <- paste(main.dir, "/RESULTS/geneExpr_vs_pval/", "geneExpr_vs_pval_JTK_", 1, ".pdf", sep="")
pdf(file = out.pdf.name, onefile=FALSE, width = 3, height = 3) # or other device
par(mfcol = c(1, 1), mar=c(0.1, 2.2, 0.1, 0.1))
smoothScatter(log2(algoTissueSubset.full.dataframe.JTK[[1]]$expr_value), -log10(algoTissueSubset.full.dataframe.JTK[[1]]$default.pvalue), xlab = NA, ylab=NA,cex = 2)
dev.off()

out.pdf.name <- paste(main.dir, "/RESULTS/geneExpr_vs_pval/", "geneExpr_vs_pval_LS_", 1, ".pdf", sep="")
pdf(file = out.pdf.name, onefile=FALSE, width = 3, height = 3) # or other device
par(mfcol = c(1, 1), mar=c(0.1, 2.2, 0.1, 0.1))
smoothScatter(log2(algoTissueSubset.full.dataframe.LS[[1]]$expr_value), -log10(algoTissueSubset.full.dataframe.LS[[1]]$default.pvalue), xlab = NA, ylab=NA,cex = 2)
dev.off()

out.pdf.name <- paste(main.dir, "/RESULTS/geneExpr_vs_pval/", "geneExpr_vs_pval_meta2d_", 1, ".pdf", sep="")
pdf(file = out.pdf.name, onefile=FALSE, width = 3, height = 3) # or other device
par(mfcol = c(1, 1), mar=c(0.1, 2.2, 0.1, 0.1))
smoothScatter(log2(algoTissueSubset.full.dataframe.meta2d[[1]]$expr_value), -log10(algoTissueSubset.full.dataframe.meta2d[[1]]$default.pvalue), xlab = NA, ylab=NA,cex = 2)
dev.off()

out.pdf.name <- paste(main.dir, "/RESULTS/geneExpr_vs_pval/", "geneExpr_vs_pval_RAIN_", 1, ".pdf", sep="")
pdf(file = out.pdf.name, onefile=FALSE, width = 3, height = 3) # or other device
par(mfcol = c(1, 1), mar=c(2.2, 2.2, 0.1, 0.1))
smoothScatter(log2(algoTissueSubset.full.dataframe.RAIN[[1]]$expr_value), -log10(algoTissueSubset.full.dataframe.RAIN[[1]]$default.pvalue), xlab = NA, ylab=NA,cex = 2)
dev.off()

for (i in 2:length(tissues.list)) {
  out.pdf.name <- paste(main.dir, "/RESULTS/geneExpr_vs_pval/", "geneExpr_vs_pval_ARS_", i, ".pdf", sep="")
  pdf(file = out.pdf.name, onefile=FALSE, width = 3, height = 3) # or other device
  par(mfcol = c(1, 1), mar=c(0.1, 0.1, 0.1, 0.1))
  smoothScatter(log2(algoTissueSubset.full.dataframe.ARS[[i]]$expr_value), -log10(algoTissueSubset.full.dataframe.ARS[[i]]$default.pvalue), xlab = NA, ylab=NA,cex = 2)
  dev.off()
  
  out.pdf.name <- paste(main.dir, "/RESULTS/geneExpr_vs_pval/", "geneExpr_vs_pval_empJTK_", i, ".pdf", sep="")
  pdf(file = out.pdf.name, onefile=FALSE, width = 3, height = 3) # or other device
  par(mfcol = c(1, 1), mar=c(0.1, 0.1, 0.1, 0.1))
  smoothScatter(log2(algoTissueSubset.full.dataframe.empJTK[[i]]$expr_value), -log10(algoTissueSubset.full.dataframe.empJTK[[i]]$default.pvalue), xlab = NA, ylab=NA,cex = 2)
  dev.off()
  
  out.pdf.name <- paste(main.dir, "/RESULTS/geneExpr_vs_pval/", "geneExpr_vs_pval_GeneCycle_", i, ".pdf", sep="")
  pdf(file = out.pdf.name, onefile=FALSE, width = 3, height = 3) # or other device
  par(mfcol = c(1, 1), mar=c(0.1, 0.1, 0.1, 0.1))
  smoothScatter(log2(algoTissueSubset.full.dataframe.GeneCycle[[i]]$expr_value), -log10(algoTissueSubset.full.dataframe.GeneCycle[[i]]$default.pvalue), xlab = NA, ylab=NA,cex = 2)
  dev.off()
  
  out.pdf.name <- paste(main.dir, "/RESULTS/geneExpr_vs_pval/", "geneExpr_vs_pval_JTK_", i, ".pdf", sep="")
  pdf(file = out.pdf.name, onefile=FALSE, width = 3, height = 3) # or other device
  par(mfcol = c(1, 1), mar=c(0.1, 0.1, 0.1, 0.1))
  smoothScatter(log2(algoTissueSubset.full.dataframe.JTK[[i]]$expr_value), -log10(algoTissueSubset.full.dataframe.JTK[[i]]$default.pvalue), xlab = NA, ylab=NA,cex = 2)
  dev.off()
  
  out.pdf.name <- paste(main.dir, "/RESULTS/geneExpr_vs_pval/", "geneExpr_vs_pval_LS_", i, ".pdf", sep="")
  pdf(file = out.pdf.name, onefile=FALSE, width = 3, height = 3) # or other device
  par(mfcol = c(1, 1), mar=c(0.1, 0.1, 0.1, 0.1))
  smoothScatter(log2(algoTissueSubset.full.dataframe.LS[[i]]$expr_value), -log10(algoTissueSubset.full.dataframe.LS[[i]]$default.pvalue), xlab = NA, ylab=NA,cex = 2)
  dev.off()
  
  out.pdf.name <- paste(main.dir, "/RESULTS/geneExpr_vs_pval/", "geneExpr_vs_pval_meta2d_", i, ".pdf", sep="")
  pdf(file = out.pdf.name, onefile=FALSE, width = 3, height = 3) # or other device
  par(mfcol = c(1, 1), mar=c(0.1, 0.1, 0.1, 0.1))
  smoothScatter(log2(algoTissueSubset.full.dataframe.meta2d[[i]]$expr_value), -log10(algoTissueSubset.full.dataframe.meta2d[[i]]$default.pvalue), xlab = NA, ylab=NA,cex = 2)
  dev.off()
  
  out.pdf.name <- paste(main.dir, "/RESULTS/geneExpr_vs_pval/", "geneExpr_vs_pval_RAIN_", i, ".pdf", sep="")
  pdf(file = out.pdf.name, onefile=FALSE, width = 3, height = 3) # or other device
  par(mfcol = c(1, 1), mar=c(2.2, 0.1, 0.1, 0.1))
  smoothScatter(log2(algoTissueSubset.full.dataframe.RAIN[[i]]$expr_value), -log10(algoTissueSubset.full.dataframe.RAIN[[i]]$default.pvalue), xlab = NA, ylab=NA,cex = 2)
  dev.off()
}


out.pdf.name <- paste(main.dir, "/RESULTS/geneExpr_vs_pval/", "geneExpr_vs_pval_RAIN_", 4, "_BIS.pdf", sep="")
pdf(file = out.pdf.name, onefile=FALSE, width = 3, height = 3) # or other device
par(mfcol = c(1, 1), mar=c(2.2, 2.2, 0.1, 0.1))
smoothScatter(log2(algoTissueSubset.full.dataframe.ARS[[4]]$expr_value), -log10(algoTissueSubset.full.dataframe.ARS[[4]]$default.pvalue), xlab = NA, ylab=NA,cex = 2)
dev.off()

smoothScatter(log2(algoTissueSubset.full.dataframe.ARS[[4]]$expr_value), -log10(algoTissueSubset.full.dataframe.ARS[[4]]$default.pvalue), postPlotHook = fudgeit)












pval.text <- gsub(".pvalue", " p-value", p.val)
#pValuesGeneExprPlot <- 
ggplot(full.dataframe, aes(x=log2(expr_value), y=-log10(full.dataframe[[p.val]]))) + 
  #geom_bin2d(bins = 130) +
  geom_point(aes(x=log2(expr_value), y=-log10(full.dataframe[[p.val]])), size=0.0001, color="gray34") +
  #stat_density2d(aes(fill = ..density..^0.2), geom = "tile") +
  scale_fill_continuous(low = "white", high = "lightgoldenrod4") +
  facet_grid(algorithm ~ tissue, scales = "free") +
  theme_light(base_line_size = 0.2) +
  labs(x = "log2(median(gene expression))", y = paste("-log10(", pval.text, ")", sep="") ) +
  theme(axis.text.x = element_text(size=9, hjust = 0.5, vjust = 0.5),
        axis.title.x = element_text(size=10, face="bold", hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(size=10, face="bold.italic", hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size=6, face="italic", hjust = 0.5, vjust = 0.5),
        strip.text.x = element_text(size = 10, face="bold", color = "brown4", hjust = 0.5, vjust = 0.5),
        strip.text.y = element_text(size = 10, face="bold", color = "darkblue", hjust = 0.5, vjust = 0.5),
        legend.position = "none") 

ggplot(data = full.dataframe, aes(x = log2(expr_value), y = -log10(default.pvalue))) +   
  stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE, n = 200) +   
  scale_fill_continuous(low = "white", high = "dodgerblue4") +
  facet_grid(algorithm ~ tissue)


#### #### #### #### #### 
#### Classic Scatter plot with a line for the smoothed conditional mean.
finalPlot <- ggplot(data = full.dataframe,  aes(x = log2(expr_value), y = -log10(default.pvalue))) +
  geom_point(shape = ".", col="grey") +
  facet_grid(algorithm ~ tissue, scales = "free_y") +
  geom_smooth(span = 0.2)

out.pdf.name <- paste(main.dir, "/RESULTS/geneExpr_vs_pval/", "tissuesComparison_", p.val, "_vs_geneExpr_", 
                      species, ".pdf", sep="")
pdf(file = out.pdf.name, onefile=FALSE, width = 5.2, height = 6.9) # or other device
print(finalPlot)
dev.off()
#### #### #### #### #### 



###  Export as pdf file
logo <- paste("~/Documents/workspace/DATA/images/", species.name, "_image.png", sep="")
out.pdf.name <- paste(main.dir, "/RESULTS/other_images/", "tissuesComparison_", p.val, "_vs_geneExpr_", 
                      species, ".pdf", sep="")
final.draw.plot <- ggdraw() +
  draw_plot(pValuesGeneExprPlot, 
            x = 0.12, y = 0, width = 0.88, height = 1) +
  draw_image(logo,
             x = 0, y = 0.82,
             width = 0.15, height = 0.15) 

width = 2+length(tissue.list)
pdf(file = out.pdf.name, onefile=FALSE, width = width, height = 6.5) # or other device
print(final.draw.plot)
dev.off()





