######################################################################
######################################################################
# log(meadian(gene expression)) VS p-values given by algorithms ######
########## Comparing tissues #######################################
######################################################################
######################################################################
library(ggplot2)
require(RColorBrewer)

main.dir <- "~/Documents/workspace"
######
species <- "anopheles"
tissue.list <- c("head_LD", "head_DD")
######

file.name <- paste(main.dir, "DATA", species, tissue.list[1], "normalized_default.pvalue/normalized_pvalue_per_gene.txt", sep = "/")
normalized.data.LD <- read.table(file.name, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
file.name <- paste(main.dir, "DATA", species, tissue.list[2], "normalized_default.pvalue/normalized_pvalue_per_gene.txt", sep = "/")
normalized.data.DD <- read.table(file.name, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (normalized.data.LD$ID == normalized.data.DD$ID) {
  merged.table <- data.frame(ID=normalized.data.DD$ID, under_LD=normalized.data.LD$pvalue_brown, under_DD=normalized.data.DD$pvalue_brown, algorithm=normalized.data.DD$algorithm)
}

# Re-order the order of algorithms
ordered.group <- c("ARS", "GeneCycle", "JTK", "LS", "RAIN", "empJTK", "meta2d")
merged.table$algorithm <- factor(merged.table$algorithm, levels=ordered.group)

# Correlation 
subset.table <- subset(merged.table, algorithm=="meta2d")
cor.test(log10(subset.table$under_LD), log10(subset.table$under_DD), na.action="na.exclude", method = "pearson")

########
# PLOT #
########
#finalPlot <- 
ggplot(merged.table, aes(x=-log10(under_LD), y=-log10(under_DD))) + 
  geom_point(aes(x=-log10(under_LD), y=-log10(under_DD)), size=0.0001, color="gray34") +
  scale_fill_continuous(low = "white", high = "lightgoldenrod4") +
  facet_wrap(~ algorithm, scales = "free", ncol = 2) +
  theme_light(base_line_size = 0.2) +
  labs(y = "-log10(p-value under LD)", x = "-log10(p-value under DD)") +
  theme(strip.text.x = element_text(color = "black", face = "bold"),
        strip.text.y = element_text(color = "black", face = "bold"))


