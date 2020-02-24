######################################################################
######################################################################
# log(meadian(gene expression)) VS p-values given by algorithms ######
########## Comparing tissues #######################################
######################################################################
######################################################################
library(tidyverse)
library(ggplot2)
require(RColorBrewer)
library(ReporteRs)
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


# Remove methods which does not get amplitude in output: 
full.dataframe <- subset(full.dataframe, algorithm %in% c("ARS", "empJTK", "JTK", "LS", "meta2d"))

########
# PLOT #
########
finalPlot <- ggplot(full.dataframe, aes(x=log2(amplitude), y=-log10(full.dataframe[[p.val]]))) + 
  geom_point(aes(x=log2(amplitude), y=-log10(full.dataframe[[p.val]])), size=0.0001, color="gray34") +
  scale_fill_continuous(low = "white", high = "lightgoldenrod4") +
  facet_grid(algorithm ~ tissue, scales = "free") +
  theme_light(base_line_size = 0.2) +
  labs(y = "-log10(default p-value)", x = "log2(amplitude)") +
  theme(strip.text.x = element_text(color = "black", face = "bold"),
        strip.text.y = element_text(color = "black", face = "bold"))






out.pdf.name <- paste(main.dir, "/RESULTS/amplitude_vs_pval.pdf", sep="")
pdf(file = out.pdf.name, onefile=FALSE, width = 4.5, height = 4.5) # or other device
print(finalPlot)
dev.off()



