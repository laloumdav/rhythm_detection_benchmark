########################################
### Density distribution of pvalues ####
## in rhythmic orthologs ###
############### VS #####################
######### the other orthologs ##############
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


args = commandArgs(trailingOnly=TRUE)
species1 <- args[1]
tissue1 <- args[2]
species2 <- args[3]
tissue2 <- args[4]
species2.algorithm <- args[5]
species2.threshold = as.numeric(args[6])
only.within.conserved.genes <- as.logical(args[7])
only.ortho.one.to.one <- as.logical(args[8])
p.val <- args[9] # should be among "raw.pvalue", "default.pvalue", or "BH.Q"




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
master$group <- "other orthologs"
master[master[[paste(species2.name, species2.algorithm, sep="_")]] <= species2.threshold
       & master[["conserved_gene"]] == TRUE, "group"] <- "rhythmic orthologs"


final.table <- master[, c(species1.algorithm.list, "group")]
final.table <- melt(final.table, id.vars = "group")
colnames(final.table) <- c("group", "algorithm", "pvalue")



etoilesFunction <- function (pvalue) {
  return(c("n.s.", "*", "**", "***", "****")[which.max(which(pvalue <= c(1, 0.05, 0.01, 0.001, 2.2e-16)))])
}


ks.pvalue.table <- data.frame()
for (i in 1:length(species1.algorithm.list)) {
  algo.values <- subset(final.table, algorithm == species1.algorithm.list[i])
  ks.test.result <- ks.test(algo.values[ algo.values$group=="rhythmic orthologs", "pvalue"], 
                       algo.values[ algo.values$group=="other orthologs", "pvalue"])
  
  ks.pvalue <- ks.test.result$p.value
  ks.Dvalue <- ks.test.result$statistic
  ks.pvalue.table[1, species1.algorithm.list[i]] <- ks.pvalue
  ks.pvalue.table[2, species1.algorithm.list[i]] <- ks.Dvalue
  
  ks.text <- etoilesFunction(ks.pvalue)
  ks.text <- paste(ks.text, paste("(D=", round(ks.Dvalue, 3), ")", sep=""), sep=" ")
  
  final.table[final.table$algorithm == species1.algorithm.list[i], "ks.text"] <- ks.text
  final.table[final.table$algorithm == species1.algorithm.list[i], "max.density"] <- max(density(algo.values$pvalue)$y)
}
signif.table <- unique(final.table[, c("algorithm", "ks.text", "max.density")])



# To conserve to same order of groups in the plot : 
ordered.group <- c("rhythmic orthologs", "other orthologs")
final.table$group <- factor(final.table$group, levels=ordered.group)
# Re-order the order of algorithms
ordered.group <- c("ARS", "GeneCycle", "JTK", "LS", "RAIN", "empJTK", "meta2d")
final.table$algorithm <- factor(final.table$algorithm, levels=ordered.group)


  
pval.text1 <- gsub(".pvalue", "\n p-value", p.val)
pval.text2 <- gsub(".pvalue", " p-value", p.val)
##########################
######### PLOT ###########
##########################
densityDistribPlot <- ggplot(final.table, aes(x = pvalue, fill = fct_rev(factor(group)), color = NA)) +
    geom_density(aes(x = pvalue, fill = factor(group)), alpha=0.7, size=0.2, color="black") +
    scale_fill_manual(values = c(`rhythmic orthologs`="blue1", `other orthologs`="sandybrown")) +
    facet_wrap(~ algorithm, scales = "free", ncol = 2) +
    geom_text(data = signif.table, 
              aes(y = 0.8*max.density, x = 0.5, label = ks.text), inherit.aes = FALSE) +
    scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), limits = c(0, 1)) +
    labs(title = "Density distribution",
         x = pval.text2, y = "density") +
    theme_ridges(grid = TRUE, line_size = 0.1) +
    theme(axis.text.x = element_text(size=9, face="bold.italic", colour="springgreen4", hjust = 0.5, vjust = 0.5),
          axis.title.x = element_text(size=10, face="bold.italic", colour="springgreen4", hjust = 0.5, vjust = 0.5),
          axis.title.y = element_text(size=10, face="bold", hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size=8, hjust = 0.5, vjust = 0.5),
          strip.text = element_text(size=12, face="bold.italic", colour="darkorchid4", hjust = 0.5, vjust = 0.5),
          strip.background = element_rect(fill="transparent", color=NA),
          plot.title = element_text(size=13, face="bold", hjust = 0.5),
          legend.position = c(0.53, 0.1),
          legend.direction = 'vertical',
          legend.title = element_blank(),
          legend.text = element_text(size=11, face="bold.italic"),
          legend.key = element_rect(color = "transparent", fill = "transparent"),
          legend.spacing.x = unit(0.3, "cm"),
          legend.spacing.y = unit(0.3, "cm"),
          axis.ticks.y = element_blank())



  
  
###  Export as pdf file
ratio.illust <- paste(main.dir, "/DATA/images/ratio_illustration.png", sep="")
logo1 <- paste(main.dir, "/DATA/images/", species1.name, "_image.png", sep="")
logo2 <- paste(main.dir, "/DATA/images/", species2.name, "_image.png", sep="")
ortho.groups.img <- paste(main.dir, "/DATA/images/ortho_groups.png", sep="")
ks.test.label <- paste("Kolmogorov-Smirnov test:\n",
                       "n.s.:  ]0.05-1]", 
                       "*:  ]0.01-0.05]",
                       "**:  ]0.001-0.01]", 
                       "***:  ]2.2e-16-0.001]", 
                       "****:  <2.2e-16", 
                       "D = Kolmogorov's D statistic", sep="\n")
out.pdf.name <- paste(main.dir, "/RESULTS/ortho_pval_distrib_images/", p.val, "_ORTHO_distrib_", 
                      species1, "_", species2, "_", tissue1, ".pdf", sep="")

final.draw.plot <- ggdraw() +
  draw_plot(densityDistribPlot, 
            x = 0.05, y = 0, width = 0.7, height = 0.6) +
  draw_label(ks.test.label, size = 9, fontfamily = "serif",
             x = 0.76, y = 0.1, hjust = 0) +
  
  draw_image(ratio.illust, scale = 1, 
             x = -0.1, y = 0.7, width = 0.7, height = 0.2) +
  draw_image(logo1, scale = 1, 
             x = 0.138, y = 0.77, width = 0.1, height = 0.1) +
  draw_image(logo2, scale = 1, 
             x = 0.25, y = 0.77, width = 0.1, height = 0.1) +
  draw_label("algorithm", colour = "darkorchid4", fontface = "bold.italic", size = 11, 
             x = 0.11, y = 0.86) +
  draw_label(pval.text1, colour = "springgreen4", fontface = "bold.italic", size = 11, 
             x = 0.115, y = 0.76) +
  draw_label(species2.algorithm, colour = "black", fontface = "bold", size = 11, 
             x = 0.39, y = 0.86) +
  draw_label(species2.threshold, colour = "black", fontface = "bold", size = 11, 
             x = 0.39, y = 0.76) +
  draw_label(tissue1, colour = "black", fontface = "bold", size = 11, 
             x = 0.21, y = 0.72) +
  draw_label(tissue2, colour = "black", fontface = "bold", size = 11, 
             x = 0.3, y = 0.72) +
  
  draw_image(ortho.groups.img, scale = 1.5, 
             x = 0.5, y = 0.7, width = 0.5, height = 0.2) +
  draw_image(logo1, scale = 1, 
             x = 0.63, y = 0.88, width = 0.11, height = 0.11) +
  draw_label(paste(species1.name, "genes", sep=" "), colour = "black", fontface = "bold", size = 12, 
             x = 0.6, y = 0.87) +
  draw_image(logo2, scale = 0.6, 
             x = 0.821, y = 0.72, width = 0.1, height = 0.1) +
  draw_label(tissue2, colour = "black", fontface = "bold", size = 9, 
             x = 0.9, y = 0.74) +
  draw_image(logo2, scale = 0.6, 
             x = 0.79, y = 0.655, width = 0.1, height = 0.1)


pdf(file = out.pdf.name, onefile=FALSE, width = 8, height = 8) # or other device
print(final.draw.plot)
dev.off()


# Also save it as a R object for R Rmakdown:
out.rds.name <- paste(main.dir, "/RESULTS/ortho_pval_distrib_images/", p.val, "_ORTHO_distrib_", 
                      species1, "_", species2, "_", tissue1, ".rds", sep="")
saveRDS(final.draw.plot, file = out.rds.name)
