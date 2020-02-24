##############################
### per gene normalization ###
##############################

#source("https://bioconductor.org/biocLite.R")
#biocLite("EmpiricalBrownsMethod")
library(EmpiricalBrownsMethod)
library(plyr)
library(metap)

main.dir <- "~/Documents/rhythm_detection_benchmark"


args = commandArgs(trailingOnly=TRUE)
species <- args[1]
tissue <- args[2]

file.dir <- paste(paste(main.dir, "DATA", species, tissue, sep = "/"), "/", sep="")

# check if working directory exist
if (dir.exists(file.dir)==FALSE){
  stop("Directory/file doesn't exist. Maybe: mouse_RNAseq or mouse_microarray ???")
}


algorithms.list <- list.files(file.dir)
algorithms.list <- algorithms.list[ grepl("empJTK|RAIN|JTK|LS|ARS|meta2d|GeneCycle", algorithms.list) ]
algorithms.list.name <- gsub(".txt", "", algorithms.list)

algorithms.list.files <- paste(file.dir, algorithms.list, sep = "")

system(paste("echo pvaluePerGeneNormalization is working with existing output, ie :", paste(algorithms.list, collapse=' , ' ), sep=" "))


file.name <- paste(file.dir, tissue, ".txt", sep = "")
raw.data <- read.table(file.name, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)

###############
### SCRIPT ####
###############
p.value <- c("raw.pvalue", "default.pvalue", "BH.Q")
for (a in 1:length(p.value)) {
  p.val <- p.value[a]

  
  algorithm <- NULL
  algorithm.normalized <- NULL
  for (i in 1:length(algorithms.list.files)){
  
    algorithm.result <- read.table(algorithms.list.files[i], head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  
    ##############################################################################################
    ##############################################################################################
    ############## if has not been previously done ##############
    ######  Ajusting all raw p-values by the same method #######
    ###### getting q-value with Benjamini-Hochberg method ######
    if ( !("BH.Q" %in% colnames(algorithm.result)) ) {
      algorithm.result$`BH.Q` <- p.adjust(algorithm.result$raw.pvalue, method = "BH")
      algorithm.result <- algorithm.result[ , c("ID", "raw.pvalue", "default.pvalue", "BH.Q", "period", "phase", "amplitude")]
      write.table(algorithm.result, algorithms.list.files[i], row.names = F, quote = F, sep = "\t")
        }
    ##############################################################################################
    ##############################################################################################
    
    
    algorithm[[i]] <- algorithm.result
    algorithm[[i]] <- algorithm[[i]][, c("ID", p.val)]
  
    # check if normalization is needed :
    if (length(unique(algorithm.result$ID)) != length(algorithm.result$ID)){
  
      # NORMALIZATION for genes having several data (ProbIDs or Transcripts) : several options : 
      # 1) Mean. Also create algorithm.normalized data.frame
      algorithm.normalized[[i]] <- aggregate(algorithm[[i]][-1], algorithm[[i]][1], mean)
      colnames(algorithm.normalized[[i]]) <- c("ID", "pvalue_mean")
  
      # 2) fisher, kost and brown normalization:
      genes.counts <- plyr::count(algorithm[[i]]$ID)
      genes.in.several.copies <- subset(genes.counts, freq >1)[,"x"]
  
      original.data <- NULL
      pvalues <- NULL
      for (k in 1:nrow(algorithm.normalized[[i]])) {
        gene.id <- algorithm.normalized[[i]]$ID[k]
  
        if (gene.id %in% genes.in.several.copies) {
          original.data[[k]] <- raw.data[grep(gene.id, raw.data$ID),-1]
          # modify the original.data[[k]] adding 1e-16 at the first column if all columns have 0 values
          # because empiricalBrownsMethod() does not work otherwise
          # => this step is not needed if these data have been removed initialy
          for (l in 1:nrow(original.data[[k]])){
            if (sum(original.data[[k]][l, ])==0){
              original.data[[k]][l, 1] <- 1e-16
            }
          }
          
          pvalues[[k]] <- algorithm[[i]][grep(gene.id, algorithm[[i]]$ID), p.val]
          empirical.browns <- empiricalBrownsMethod(data_matrix=original.data[[k]], p_values=pvalues[[k]], extra_info=TRUE)
          empirical.kosts <- kostsMethod(data_matrix=as.matrix(original.data[[k]]), p_values=pvalues[[k]], extra_info=TRUE)
          metap <- allmetap(pvalues[[k]], method = "all")
          
          pvalue.brown <- empirical.browns[[1]]
          pvalue.fisher <- empirical.browns[[2]]
          pvalue.kost <- empirical.kosts[[1]]
          pvalue.logitp <- unlist(metap["logitp", "p"], use.names = FALSE)
          pvalue.sumlog <- unlist(metap["sumlog", "p"], use.names = FALSE)
          pvalue.sump <- unlist(metap["sump", "p"], use.names = FALSE)
          pvalue.sumz <- unlist(metap["sumz", "p"], use.names = FALSE)
  
          algorithm.normalized[[i]]$pvalue_brown[k] <- pvalue.brown
          algorithm.normalized[[i]]$pvalue_fisher[k] <- pvalue.fisher
          algorithm.normalized[[i]]$pvalue_kost[k] <- pvalue.kost
          algorithm.normalized[[i]]$pvalue_logitp[k] <- pvalue.logitp
          algorithm.normalized[[i]]$pvalue_sumlog[k] <- pvalue.sumlog
          algorithm.normalized[[i]]$pvalue_sump[k] <- pvalue.sump
          algorithm.normalized[[i]]$pvalue_sumz[k] <- pvalue.sumz
  
        } else {
          algorithm.normalized[[i]]$pvalue_brown[k] <- algorithm[[i]][grep(gene.id, algorithm[[i]]$ID), p.val]
          algorithm.normalized[[i]]$pvalue_fisher[k] <- algorithm[[i]][grep(gene.id, algorithm[[i]]$ID), p.val]
          algorithm.normalized[[i]]$pvalue_kost[k] <- algorithm[[i]][grep(gene.id, algorithm[[i]]$ID), p.val]
          algorithm.normalized[[i]]$pvalue_logitp[k] <- algorithm[[i]][grep(gene.id, algorithm[[i]]$ID), p.val]
          algorithm.normalized[[i]]$pvalue_sumlog[k] <- algorithm[[i]][grep(gene.id, algorithm[[i]]$ID), p.val]
          algorithm.normalized[[i]]$pvalue_sump[k] <- algorithm[[i]][grep(gene.id, algorithm[[i]]$ID), p.val]
          algorithm.normalized[[i]]$pvalue_sumz[k] <- algorithm[[i]][grep(gene.id, algorithm[[i]]$ID), p.val]
        }
      }
  
    } else {
      system(paste("echo Normalization not needed. Each gene has a unique data. A final file is created anyway, called pvalue_per_gene.txt"))
      algorithm.normalized[[i]] <- algorithm[[i]]
    }
  
    algorithm.normalized[[i]]$algorithm <- rep(algorithms.list.name[i], nrow(algorithm.normalized[[i]]))
  }
  
  
  # We make a unique dataframe with these data:
  full.dataframe <- data.frame()
  for (i in 1:length(algorithms.list)){
    full.dataframe <- rbind(full.dataframe, algorithm.normalized[[i]])
  }
  
  if (length(unique(algorithm.result$ID)) != length(algorithm.result$ID)){
    colnames(full.dataframe) <- c("ID", "pvalue_mean", "pvalue_brown", "pvalue_fisher", "pvalue_kost", "pvalue_logitp", "pvalue_sumlog", "pvalue_sump" , "pvalue_sumz", "algorithm")
  } else {
    colnames(full.dataframe) <- c("ID", "pvalue", "algorithm")
  }
  
  
  normalized.pval.dir <- paste(file.dir, "normalized_", p.val, sep="")
  if (!dir.exists(normalized.pval.dir)){
    dir.create(file.path(normalized.pval.dir))
  }
  write.table(full.dataframe, paste(normalized.pval.dir, "normalized_pvalue_per_gene.txt", sep="/"), row.names = F, quote = F, sep = "\t")

}
