######################################################################
######################################################################
### Gene expression variation over time, of some genes ###############
######################################################################
######################################################################
library(ggplot2)
require(RColorBrewer)
library(esquisse)
library(ggplot2)

main.dir <- "~/Documents/rhythm_detection_benchmark"




##############################
species <- "rat"
tissue <- "lung"
##############################
if (grepl("mouse", species)){ species.name <- "mouse" } else { species.name <- species }



tissue.file <- paste(tissue, ".txt", sep="")

file.dir <- paste(paste(main.dir, "DATA", species, tissue, sep = "/"), "/", sep="")
file.name <- paste(file.dir, tissue.file, sep = "")
raw.dataset <- read.table(file.name, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)


############################################################
############################################################
if ( species.name == "mouse" | species.name == "drosophila" ) {
  circad.entrainment.genes <- read.table(paste(main.dir, "DATA", species, "kegg_circadian_rhythmc_genes.txt", sep="/"), h=T)
} else {
  ## Check if known_circadian_genes list exist 
  circad.genes.file <- paste(main.dir, "DATA", species, "known_circadian_genes.txt", sep="/")
  if (file.exists(circad.genes.file)) {
    circad.entrainment.genes <- read.table(circad.genes.file, h=T)
  } else { stop("!! Circadian genes list does not exist !!") }
}

if (ncol(circad.entrainment.genes) <2) {
  circad.entrainment.genes$ID <- circad.entrainment.genes[,1]
  colnames(circad.entrainment.genes)[1] <- "Gene.Name"
}

colnames(circad.entrainment.genes)[2] <- "ID"

ID.circad.entrainment.genes <- as.data.frame(circad.entrainment.genes$ID)
colnames(ID.circad.entrainment.genes) <- "ID"


# keep only these genes : 
raw.dataset <- merge(ID.circad.entrainment.genes, raw.dataset, by="ID")
# Remove duplicated data for one same gene :
raw.dataset <- raw.dataset[duplicated(raw.dataset$ID) == FALSE, ]
# Retreive Gene.Name
raw.dataset <- merge(circad.entrainment.genes, raw.dataset, by="ID")
raw.dataset <- raw.dataset[, -1]

restricted.raw.dataset <- raw.dataset

# stack function manner : 
final.table <- stack(restricted.raw.dataset, select = -Gene.Name)
final.table$Gene.Name <- rep(restricted.raw.dataset$Gene.Name, ncol(restricted.raw.dataset)-1)
colnames(final.table) <- c("expr_value", "time_point", "Gene.Name")
final.table$time_point <- gsub("CT|ZT", "", final.table$time_point)
final.table$time_point <- as.numeric(final.table$time_point)


#nb.col=3

# Plot
ggplot.plot <- 
  ggplot(final.table, aes(x = time_point, y = expr_value)) +
  geom_line(color = '#0c4c8a') +
  geom_smooth(aes(x=time_point, y=expr_value, col=Gene.Name), color = "#FC4E07", 
              span = 0.5, se = FALSE, size = 0.7, method = "loess") +
  #theme_linedraw() +
  labs(x = 'time-point (hours)', y = 'gene expression') +
  facet_wrap( ~ Gene.Name, scales = "free") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y =element_text(size=6),
        axis.title.x = element_text(face="bold", hjust=0.5),
        axis.title.y = element_text(face="bold", hjust=0.5),
        legend.text = element_text(size = 11),
        legend.background = element_rect(color = "black", linetype = "solid", size = 0.2),
        legend.key = element_rect(colour = "white", fill = NA), 
        legend.title = element_text(face="bold", hjust=0.5, size=12),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

ggplot.plot.2 <- 
  ggplot(final.table, aes(x = time_point, y = expr_value)) +
  geom_line(color = '#0c4c8a') +
  labs(x = 'time-point (hours)', y = 'gene expression') +
  facet_wrap( ~ Gene.Name, scales = "free") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y =element_text(size=6),
        axis.title.x = element_text(face="bold", hjust=0.5),
        axis.title.y = element_text(face="bold", hjust=0.5),
        legend.text = element_text(size = 11),
        legend.background = element_rect(color = "black", linetype = "solid", size = 0.2),
        legend.key = element_rect(colour = "white", fill = NA), 
        legend.title = element_text(face="bold", hjust=0.5, size=12),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
  

###  Export as pdf file
logo <- paste(main.dir, "/DATA/images/", species.name, "_image.png", sep="")
tissue.text <- gsub("_", " ", tissue)

### 
# Plot 1
###
out.pdf.name <- paste(main.dir, "/RESULTS/other_images/", "geneExpr_over_time", 
                      species, "_", tissue, "_fittingCurve.pdf", sep="")
final.draw.plot <- ggdraw() +
  draw_plot(ggplot.plot, 
            x = 0.12, y = 0, width = 0.88, height = 1) +
  draw_image(logo,
             x = 0, y = 0.82,
             width = 0.15, height = 0.15) +
  draw_label(tissue.text, colour = "black", fontface = "bold", size = 10, 
             x = 0.09, y = 0.81)

width = 0.3*nrow(circad.entrainment.genes)
pdf(file = out.pdf.name, onefile=FALSE, width = 8.5, height = 6)
print(final.draw.plot)
dev.off()

### 
# Plot 2
###
out.pdf.name <- paste(main.dir, "/RESULTS/other_images/", "geneExpr_over_time", 
                      species, "_", tissue, ".pdf", sep="")
final.draw.plot <- ggdraw() +
  draw_plot(ggplot.plot.2, 
            x = 0.12, y = 0, width = 0.88, height = 1) +
  draw_image(logo,
             x = 0, y = 0.82,
             width = 0.15, height = 0.15) +
  draw_label(tissue.text, colour = "black", fontface = "bold", size = 10, 
             x = 0.09, y = 0.81)

#width = 0.4*nrow(circad.entrainment.genes)
pdf(file = out.pdf.name, onefile=FALSE, width = 8.5, height = 6)
print(final.draw.plot)
dev.off()

