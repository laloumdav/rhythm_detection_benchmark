########################################
### Density distribution of pvalues ####
## in conserved rhythmic genes group ###
############### VS #####################
######### the other genes ##############
########################################
library(ggplot2)
library(reshape2)
require(RColorBrewer)
library(MASS)
library(aggregation)
library(magick)
library(ggridges)
library(cowplot)

main.dir <- "~/Documents/rhythm_detection_benchmark"



species1 <- "mouse_microarray"
tissue1 <- "lung"
species2 <- "rat"
tissue2 <- "lung"
species2.algorithm <- "ARS"
#species2.threshold = as.numeric(args[6])
only.within.conserved.genes <- TRUE
only.ortho.one.to.one <- TRUE
p.val <- "default.pvalue"




if (grepl("mouse", species1)){ species1.name <- "mouse" } else { species1.name <- species1 }
if (grepl("mouse", species2)){ species2.name <- "mouse" } else { species2.name <- species2 }


file.name <- paste(p.val, species1, tissue1, species2, tissue2, sep="_")
orthology.per.gene <- read.table(paste(main.dir, "/DATA/ORTHOLOGY/result/", file.name, ".txt", sep=""), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# Retrieve pvalues of species1 for each algorithm :
species1.pvalues <- dcast(orthology.per.gene,  
                          orthology.per.gene[, paste(species1.name, "ID", sep="_")] ~ orthology.per.gene[, paste(species1.name, "algorithm", sep="_")], 
                          value.var = paste(species1.name, "pvalue", sep="_"), 
                          fun.aggregate = mean, drop=TRUE, fill=1)
colnames(species1.pvalues)[1] <- paste(species1.name, "ID", sep="_")
# We replace NA by value=1 because it correspond to outsider genes with absolutely NO signal of rhythmicity
species1.pvalues <- replace(species1.pvalues, is.na(species1.pvalues), 1)

# List of algorithms for which there is data :
species1.algorithm.list <- colnames(species1.pvalues)[colnames(species1.pvalues) != paste(species1.name, "ID", sep="_")]




### CONSERVATION info ###
conservation <- unique(orthology.per.gene[ , c(paste(species1.name, "ID", sep="_"), "orthotype")])
#### Orthology ####
orthologs.pvalues <- dcast(orthology.per.gene,  
                           orthology.per.gene[, paste(species1.name, "ID", sep="_")] ~ orthology.per.gene[, paste(species2.name, "algorithm", sep="_")],
                           value.var = paste(species2.name, "pvalue", sep="_"), 
                           fun.aggregate = mean, drop=TRUE, fill=NULL)
# Remove NA column if exist
orthologs.pvalues <- orthologs.pvalues[,colnames(orthologs.pvalues)!="NA"]
column.names <- colnames(orthologs.pvalues)
column.names <- c(paste(species1.name, "ID", sep="_"), paste(species2.name, column.names[-1], sep="_"))
colnames(orthologs.pvalues) <- column.names
# Is there data for the species2 ortholog gene :
orthologs.pvalues$species2_data_available <- apply(orthologs.pvalues[,-1], 1, function(x) all(!is.nan(x)))
# Conserved genes: Affect TRUE if the gene of species1 has one or more orthologs in species2
orthologs.pvalues <- merge(orthologs.pvalues, conservation, by=paste(species1.name, "ID", sep="_"))
orthologs.pvalues$conserved_gene <- !is.na(orthologs.pvalues$orthotype)




# Remove conserved genes that do not have data in species2 :
master <- subset(orthologs.pvalues, !(species2_data_available==FALSE & conserved_gene==TRUE))

# If argument only.within.conserved.genes = TRUE : we keep only conserved genes to work with
if (isTRUE(only.within.conserved.genes)){
  master <- subset(master, conserved_gene == TRUE)
}
# If argument only.ortho.one.to.one = TRUE : we keep only orthologs 1:1
if (isTRUE(only.ortho.one.to.one)){
  master <- subset(master, orthotype == "1:1")
}
master <- merge(master, species1.pvalues, by.x=paste(species1.name, "ID", sep="_"))


ordered.group <- c("ARS", "GeneCycle", "JTK", "LS", "RAIN", "empJTK", "meta2d")
final.table <- data.frame()
for (i in 1: length(ordered.group)) {
  algo.name <- ordered.group[i]
  algo.name.species2 <- paste(species2.name, algo.name, sep="_")
  species.ID <- paste(species1.name, "ID", sep="_")
  subset.master <- master[, c(species.ID, algo.name, algo.name.species2)]
  colnames(subset.master) <- c(species.ID, "pvalue_species1", "pvalue_species2")
  subset.master$algorithm <- rep(paste(algo.name, algo.name, sep=" - "), nrow(subset.master))
  
  final.table <- rbind(final.table, subset.master)
}

# Re-order the order of algorithms
ordered.group <- paste(ordered.group, ordered.group, sep=" - ")
final.table$algorithm <- factor(final.table$algorithm, levels=ordered.group)


finalPlot <- ggplot(final.table, aes(x=-log10(pvalue_species1), y=-log10(pvalue_species2))) + 
  geom_point(aes(x=-log10(pvalue_species1), y=-log10(pvalue_species2)), size=0.0001, color="gray34") +
  facet_wrap(~ algorithm, scales = "free", ncol = 2) +
  theme_light(base_line_size = 0.2) +
  theme_ridges(font_size = 12, grid = TRUE, line_size = 0.1) +
  labs(x = "-log10(mouse default p-value)", y = "-log10(rat default p-value)") +
  theme(axis.text.x = element_text(size=9, hjust = 0.5, vjust = 0.5),
      axis.text.y = element_text(size=9, hjust = 0.5, vjust = 0.5),
      axis.title.x = element_text(size=10, face="italic", hjust = 0.5, vjust = 0.5),
      axis.title.y = element_text(size=10, face="italic", hjust = 0.5, vjust = 0.5),
      strip.text.x = element_text( size=11, face="bold", hjust = 0.5, vjust = 0.5),
      plot.title = element_blank()) 
#
###  Export as pdf file
out.pdf.name <- paste(main.dir, "/RESULTS/", "orthologsPval_mouse_vs_rat.pdf", sep="")
pdf(file = out.pdf.name, onefile=FALSE, width = 5, height = 8) # or other device
print(finalPlot)
dev.off()



subset.final.table <- subset(final.table, algorithm == "LS - LS")
cor.test(asinh(subset.final.table$pvalue_species1), asinh(subset.final.table$pvalue_species2), na.action="na.exclude", method = "pearson")

