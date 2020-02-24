######################################
### Nb of rhythmic genes detected ####
######################################
#install.packages("venn")
library(venn)
library(UpSetR)
library(wesanderson)
library(plyr)
library(ggplot2)
library(reshape2)

main.dir <- "~/Documents/rhythm_detection_benchmark"

args = commandArgs(trailingOnly=TRUE)
species <- args[1]
tissue <- args[2]
threshold.1 <- as.numeric(args[3])
threshold.2 <- as.numeric(args[4])
p.adj.method <- args[5]


if (threshold.2 > threshold.1) {
  min.threshold <- threshold.1
  max.threshold <- threshold.2
  threshold.1 <- max.threshold
  threshold.2 <- min.threshold
}

if (grepl("mouse", species)){ species.name <- "mouse" } else { species.name <- species }

# normalization methods : PARAMETERS 
selected.methods.for.normalization <- read.csv(paste(main.dir, "scripts/normalization_parameters.csv", sep="/"), nrows = 2, sep = ",") # which type of normalization for genes with several data (microarray ProbIDs or TranscriptIDs)
per.gene.normalization.method.name <- as.character(selected.methods.for.normalization[selected.methods.for.normalization$normalization_step == "per_gene", "normalization_method"])
per.gene.normalization.method <- paste("pvalue", per.gene.normalization.method.name, sep="_")


file.dir <- paste(paste(main.dir, "DATA", species, tissue, sep = "/"), "/", sep="")

pvalues <- c("raw.pvalue", "default.pvalue", "BH.Q")
for (a in 1:length(pvalues)) {
  p.val <- pvalues[a]
  file.name <- paste(file.dir, "normalized_", p.val, "/normalized_pvalue_per_gene.txt", sep = "")
  
  normalized.file <- read.table(file.name, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  # If initial normalization per gene was needed : we keep only the one selected by pvalue.per.gene.normalization.method argument
  if (length(colnames(normalized.file))>3){
    normalized.file <- normalized.file[,c("ID", per.gene.normalization.method, "algorithm")]
  }
  colnames(normalized.file) <- c("ID", "pvalue", "algorithm")
  
  
  tot.genes.nb = length(unique(normalized.file$ID))
  list.genes.per.algorithm <- split(normalized.file, f = normalized.file$algorithm)
  
  # keep only rhythmic genes list 
  # with p-val < threshold.1
  list.rhythmicgenes.per.algorithm_threshold.1 <- list()
  for (i in 1:length(list.genes.per.algorithm)){
    data <- list.genes.per.algorithm[[i]]
    
    algorithm.name <- data[1, "algorithm"]
    # adjust pvalues according to p.adj.method
    # can be : "none" , "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"
    data$pvalue <- p.adjust(data$pvalue, method = p.adj.method)
    # subset of rhythmic genes with cutoff < threshold
    data <- subset(data, pvalue <= threshold.1)
    list.rhythmicgenes.per.algorithm_threshold.1[[algorithm.name]] <- unique(data$ID)
  }
  # with p-val < threshold.2
  list.rhythmicgenes.per.algorithm_threshold.2 <- list()
  for (i in 1:length(list.genes.per.algorithm)){
    data <- list.genes.per.algorithm[[i]]
    
    algorithm.name <- data[1, "algorithm"]
    # adjust pvalues according to p.adj.method
    # can be : "none" , "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"
    data$pvalue <- p.adjust(data$pvalue, method = p.adj.method)
    # subset of rhythmic genes with cutoff < threshold
    data <- subset(data, pvalue <= threshold.2)
    list.rhythmicgenes.per.algorithm_threshold.2[[algorithm.name]] <- unique(data$ID)
  }
  
  # Take genes kept with the highest threshold
  genes <- unique(unlist(list.rhythmicgenes.per.algorithm_threshold.1))
  
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
  
  intersect.table.rhythmicgenes_threshold.1 <- FromList(list.rhythmicgenes.per.algorithm_threshold.1, genes)
  intersect.table.rhythmicgenes_threshold.2 <- FromList(list.rhythmicgenes.per.algorithm_threshold.2, genes)
  
  # Is a vector similar to another one ? 
  # The function affect "yes" or "no" for each row
  isSimilarRows <- function(table.1, table.2) {
    if (nrow(table.1) != nrow(table.2)) { stop("table.1 and table.2 must have the same number of rows") }
    
    similarity <- all(table.1[1, ] == table.2[1, ])
    for (j in 2:nrow(table.1)) {
      similarity.row <- all(table.1[j, ] == table.2[j, ])
      similarity <- c(similarity, similarity.row)
    }
    similarity[similarity == TRUE] <- "yes"
    similarity[similarity == FALSE] <- "no"
    return(similarity)
  }
  
  
  # Is the combination conserved for the small p.value threshold ?
  algo.combinations <- isSimilarRows(intersect.table.rhythmicgenes_threshold.1, intersect.table.rhythmicgenes_threshold.2)
  intersect.table.rhythmicgenes <- intersect.table.rhythmicgenes_threshold.1
  intersect.table.rhythmicgenes$comb.retrieved.smaller.threshold <- algo.combinations

  
  
  ################
  #  Upset plot  #
  ################
  colors.list <- c(wes_palette("Darjeeling1"), wes_palette("GrandBudapest2"))
  
  # function to draw colored histogram corresponding with the smaller p-value threshold
  upsetFunc <- function(input.data) {
    data <- (input.data["comb.retrieved.smaller.threshold"] == "yes")
  }
  # metadata to get colored line for each algorithm
  metadata <- data.frame(sets=names(list.genes.per.algorithm), algorithm=names(list.genes.per.algorithm))
  
  # Queries for upset 
  upset.query.1 <- list(query = upsetFunc, 
                        color = "seashell4", active = T, 
                        query.name = paste("cutoff=", threshold.2, sep=""))
  
  upset.query.2 <- list(query = intersects, 
                        params = list("LS", "JTK", "empJTK", "RAIN", "GeneCycle", "meta2d", "ARS"), 
                        color = "cyan4", active = F, 
                        query.name = "intersection with all algorithms")
  
  simple.intersect.table.rhythmicgenes <- intersect.table.rhythmicgenes[, colnames(intersect.table.rhythmicgenes) %in% c("LS", "JTK", "empJTK", "RAIN", "GeneCycle", "meta2d", "ARS")]
  simple.intersect.table.rhythmicgenes$algo.nb <- apply(simple.intersect.table.rhythmicgenes, 1, sum)
  count.algo.nb <- plyr::count(simple.intersect.table.rhythmicgenes$algo.nb)
  
  if (nrow(count.algo.nb) < length(list.genes.per.algorithm)) {
    upset.queries.list <- list(upset.query.1)
  } else { upset.queries.list <- list(upset.query.1, upset.query.2) }
  
  
  
  #######
  #######
  # Jaccard index heatmap
  Jaccard = function (x, y) {
    M.11 = sum(x == 1 & y == 1)
    M.10 = sum(x == 1 & y == 0)
    M.01 = sum(x == 0 & y == 1)
    return (M.11 / (M.11 + M.10 + M.01))
  }
  
  input.variables = intersect.table.rhythmicgenes[,c(-8,-9)]

  jaccard.index.table = matrix(data = NA, nrow = length(input.variables), ncol = length(input.variables))
  for (r in 1:length(input.variables)) {
    for (c in 1:length(input.variables)) {
      if (c == r) {
        jaccard.index.table[r,c] = 1
      } else if (c > r) {
        jaccard.index.table[r,c] = Jaccard(input.variables[,r], input.variables[,c])
      }
    }
  }
  
  variable.names = colnames(intersect.table.rhythmicgenes)[1:7]
  jaccard.index.table <- as.data.frame(jaccard.index.table)
  colnames(jaccard.index.table) = variable.names
  rownames(jaccard.index.table) = variable.names  
  jaccard.index.table$algorithm <- row.names(jaccard.index.table)
  
  melted_jaccard.index.table <- melt(jaccard.index.table)
  melted_jaccard.index.table[is.na(melted_jaccard.index.table)] <- 0
  
  jaccardHeatmap <- ggplot(data = melted_jaccard.index.table, aes(x=variable, y=algorithm, fill=value)) + 
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "white", high = "deepskyblue4", mid = "deepskyblue",
                         midpoint = 0.5, limit = c(0,1), space = "Lab",
                         name="Jaccard\nindex") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_text(size=14, face="bold", hjust = 0.5, vjust = 0.5)
          )
  
  #######
  #######

  
  

  pfd.name <- paste("upset_DEGREE_", p.val, "_", species, "_", tissue, ".pdf", sep="")
  pfd.name <- paste(main.dir, "/RESULTS/UpSet_images/", pfd.name, sep="")
  pdf(file = pfd.name, onefile=FALSE) # or other device
  print(
    upset(intersect.table.rhythmicgenes, 
        query.legend = "top",
        mainbar.y.label = paste("Intersection size \n cutoff = ", threshold.1, sep=""),
        nsets=length(list.genes.per.algorithm), 
        nintersects=60,
        order.by = c("freq", "degree"), decreasing = c(TRUE,TRUE),
        set.metadata = list(data = metadata, plots = list(list(type = "matrix_rows", column = "algorithm", colors = c(JTK=colors.list[1], LS=colors.list[6], empJTK=colors.list[7], ARS=colors.list[4], RAIN=colors.list[5], meta2d=colors.list[2], GeneCycle="blue"), alpha = 0.5))), 
        queries = upset.queries.list)
    )
  dev.off()
  
  pfd.name <- paste("upset_FREQ_", p.val, "_", species, "_", tissue, ".pdf", sep="")
  pfd.name <- paste(main.dir, "/RESULTS/UpSet_images/", pfd.name, sep="")
  pdf(file = pfd.name, onefile=FALSE) # or other device
  print(
    upset(intersect.table.rhythmicgenes, 
        query.legend = "top",
        mainbar.y.label = paste("Intersection size \n cutoff = ", threshold.1, sep=""),
        nsets=length(list.genes.per.algorithm), 
        nintersects=60,
        order.by = "freq",
        set.metadata = list(data = metadata, plots = list(list(type = "matrix_rows", column = "algorithm", colors = c(JTK=colors.list[1], LS=colors.list[6], empJTK=colors.list[7], ARS=colors.list[4], RAIN=colors.list[5], meta2d=colors.list[2], GeneCycle="blue"), alpha = 0.5))), 
        queries = upset.queries.list)
    )
  dev.off()

  
  ################
  # Venn diagram #
  ################
  
  #title <- paste(species, tissue, sep = " ")
  #subtitle.1 <- paste(p.val, "<", threshold.1, sep = " ")
  #subtitle.2 <- paste(p.val, "<", threshold.2, sep = " ")
  #tot.genes <- paste("tot genes =", tot.genes.nb, sep = " ")
  #p.adj <- paste("p.adj.method =", p.adj.method, sep = " ")
  
  ## order for colors for the venn diagram
  #col.algo <- colnames(intersect.table.rhythmicgenes_threshold.1)
  #col.algo <- col.algo[col.algo %in% c("ARS", "empJTK", "GeneCycle", "JTK", "LS","meta2d", "RAIN")]
  #col.algo[col.algo == "ARS"] <- colors.list[4]
  #col.algo[col.algo == "JTK"] <- colors.list[1]
  #col.algo[col.algo == "LS"] <- colors.list[6]
  #col.algo[col.algo == "empJTK"] <- colors.list[7]
  #col.algo[col.algo == "RAIN"] <- colors.list[5]
  #col.algo[col.algo == "meta2d"] <- colors.list[2]
  #col.algo[col.algo == "GeneCycle"] <- "blue"
  

  ###  Exit to powerpoint file
  #Plot_Documents.pptx <- pptx(template = paste(main.dir, '/RESULTS/UpSet_images/UpSet_plots.pptx', sep=""))
  #Plot_Documents.pptx <- addSlide(Plot_Documents.pptx, "Two Content" )
  #Plot_Documents.pptx <- addPlot(Plot_Documents.pptx, fontname_sans = "Helvetica", 
  #                               function() print( 
  #                                 venn(intersect.table.rhythmicgenes_threshold.1[colnames(intersect.table.rhythmicgenes_threshold.1) != "elements"], 
  #                                      ilabels = TRUE, 
  #                                      zcolor = col.algo, 
  #                                      opacity = 0.5,
  #                                      cexsn = 0.9, 
  #                                      cexil=1)))
  #Plot_Documents.pptx <- addPlot(Plot_Documents.pptx, fontname_sans = "Helvetica", 
  #                               function() print( 
  #                                 venn(intersect.table.rhythmicgenes_threshold.2[colnames(intersect.table.rhythmicgenes_threshold.2) != "elements"], 
  #                                      ilabels = TRUE, 
  #                                      zcolor = col.algo, 
  #                                      opacity = 0.5,
  #                                      cexsn = 0.9, 
  #                                      cexil=1)))
  #Plot_Documents.pptx = addTitle( Plot_Documents.pptx, paste(title, subtitle.1, "then", subtitle.2, tot.genes, p.adj, sep = " / "))
  #writeDoc(Plot_Documents.pptx,  paste(main.dir, '/RESULTS/UpSet_images/UpSet_plots.pptx', sep=""))

}
  

