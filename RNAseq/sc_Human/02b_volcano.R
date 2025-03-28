#===============================================================================
# Load packages, set wd and load data
#===============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(viridis)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)


# Set work directory
setwd("~/Documents/Collaborations/Brickman/61224/11_niakan_adata_combined")

# Load all_res
all_res <- read.csv("EPI_Ulix_vs_DMSO_allDEGs.csv")

#===============================================================================
# Volcano plot of all DEGs
#===============================================================================

#===============================================================================
# Volcano plot EPI DEGs with top 10 HVG
#===============================================================================


padj_threshold <-  0.05
lfc_threshold <- 1.5

# Add significance categories
all_res$sig_category <- "Not Significant"
all_res$sig_category[all_res$padj < padj_threshold & all_res$log2FoldChange > lfc_threshold] <- "Upregulated"
all_res$sig_category[all_res$padj < padj_threshold & all_res$log2FoldChange < -lfc_threshold] <- "Downregulated"

# Calculate gene counts
total_genes <- nrow(all_res)
up_genes <- sum(all_res$sig_category == "Upregulated", na.rm = TRUE)
down_genes <- sum(all_res$sig_category == "Downregulated", na.rm = TRUE)
sig_genes <- up_genes + down_genes

# Create gene count labels
stat_label <-sprintf("Total: %d\nUp: %d\nDown: %d\nPadj < 0.05 and log2FC > 1.5", 
                      total_genes, up_genes, down_genes)


sig_genes <- subset(all_res, padj < 0.05 & abs(log2FoldChange) > 1.5)
top_genes <- sig_genes[order(sig_genes$padj)[1:min(50, nrow(sig_genes))], ]
top_genes <- as.list(top_genes$X)

p2 <- ggplot(all_res, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(colour = "grey", alpha = 0.1) + # Colour all points grey
  geom_point(data = all_res[all_res$log2FoldChange < -1.5 & all_res$padj < 0.05,], 
             aes(x = log2FoldChange, y = -log10(padj)), colour = "purple", alpha = 0.1) + # Colour downreg, significant in blue
  geom_point(data = all_res[all_res$log2FoldChange > 1.5 & all_res$padj < 0.05,], 
             aes(x = log2FoldChange, y = -log10(padj)), colour = "orange", alpha = 0.1) + # Colour upreg, significant in red
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black", alpha = 0.5) + 
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "black",  alpha = 0.5) + 
  geom_text_repel(data = all_res[all_res$gene_symbol %in% top_genes & abs(all_res$log2FoldChange) > 0,], 
                  aes(x = log2FoldChange, y = -log10(padj), label = gene_symbol),
                  max.overlaps = Inf,  # Allow all labels to be displayed
                  box.padding = 0.1,     # Increase spacing between labels
                  point.padding = 0.1,  # Increase spacing from points
                  min.segment.length = 0.1,  # Minimum line segment length
                  force = 0.1,          # Increase force of repulsion
                  seed = 123,          # Set seed for reproducibility
                  size = 2,            # Smaller text size
                  segment.size = 0.1,  # Thinner connection lines
                  segment.alpha = 0.8) +  # Slightly transparent line
  labs(title = "Epiblast: ERKi vs Control", subtitle = stat_label)+
   theme_classic() + theme(text = element_text(size = 8), aspect.ratio = 0.75)
plot(p2)
ggsave("volcano_top100padj.pdf", p2)


sig_genes <- subset(all_res, pvalue < 0.05 & abs(log2FoldChange) > 1.5)
top_genes <- sig_genes[order(sig_genes$padj)[1:min(100, nrow(sig_genes))], ]
top_genes <- as.list(top_genes$X)

p2 <- ggplot(all_res, aes(x = log2FoldChange, y = -log10(pvalue))) + 
  geom_point(colour = "grey", alpha = 0.1) + # Colour all points grey
  geom_point(data = all_res[all_res$log2FoldChange < 0 & all_res$pvalue <= 0.05,], 
             aes(x = log2FoldChange, y = -log10(pvalue)), colour = "purple", alpha = 0.1) + # Colour downreg, significant in blue
  geom_point(data = all_res[all_res$log2FoldChange > 0 & all_res$pvalue <= 0.05,], 
             aes(x = log2FoldChange, y = -log10(pvalue)), colour = "orange", alpha = 0.1) + # Colour upreg, significant in red
  geom_text_repel(data = all_res[all_res$gene_symbol %in% top_genes & abs(all_res$log2FoldChange) > 0,], 
                  aes(x = log2FoldChange, y = -log10(pvalue), label = gene_symbol),
                  max.overlaps = Inf,  # Allow all labels to be displayed
                  box.padding = 0.1,     # Increase spacing between labels
                  point.padding = 0.1,  # Increase spacing from points
                  min.segment.length = 0.1,  # Minimum line segment length
                  force = 0.1,          # Increase force of repulsion
                  seed = 123,          # Set seed for reproducibility
                  size = 2,            # Smaller text size
                  segment.size = 0.1,  # Thinner connection lines
                  segment.alpha = 0.8) +  # Slightly transparent line
  labs(title = "Ulixertinib vs DMSO", subtitle = "DESeq2 pvalue < 0.05")+
  theme_classic() + theme(aspect.ratio = 0.75)

plot(p2)
ggsave("volcano_top100pvalue.pdf", p2)



#===============================================================================
### ICM
#===============================================================================

# Load ICM all_res
all_res <- read.csv("ICM_Ulix_vs_DMSO_allDEGs.csv")

#===============================================================================
# Volcano plot of all DEGs
#===============================================================================

#===============================================================================
# Volcano plot EPI DEGs with genes of interest
#===============================================================================

ICM_GOI <- read.csv("ICM_genelist.csv", header = FALSE)
genes <- as.list(ICM_GOI$V1)


# Gene list from Laura, few are sig. in most recent batch

# genes <- c("DNMT3L", "DUSP6", "ETV5")



p2 <- ggplot(all_res, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(colour = "grey", alpha = 0.1) + # Colour all points grey
  geom_point(data = all_res[all_res$log2FoldChange < 0 & all_res$padj <= 0.05,], 
             aes(x = log2FoldChange, y = -log10(padj)), colour = "purple", alpha = 0.1) + # Colour downreg, significant in blue
  geom_point(data = all_res[all_res$log2FoldChange > 0 & all_res$padj <= 0.05,], 
             aes(x = log2FoldChange, y = -log10(padj)), colour = "orange", alpha = 0.1) + # Colour upreg, significant in red
  geom_text_repel(data = all_res[all_res$gene_symbol %in% genes & all_res$log2FoldChange >0,], 
                  aes(x = log2FoldChange, y = -log10(padj), label = gene_symbol), 
                  min.segment.length = unit(0, 'lines'), nudge_y = 1, nudge_x = 2, size = 3, fontface = "italic", direction = "y") +
  geom_text_repel(data = all_res[all_res$gene_symbol %in% genes & all_res$log2FoldChange <0,], 
                  aes(x = log2FoldChange, y = -log10(padj), label = gene_symbol), 
                  min.segment.length = unit(0, 'lines'), nudge_y = 1, nudge_x = -2, size = 3, fontface = "italic", direction = "y") +
  labs(title = "Ulixertinib vs DMSO", subtitle = "DESeq2 padj < 0.05")+
  theme_classic() + theme(aspect.ratio = 1)

plot(p2)
ggsave("ICM_volcano_GOI.pdf", p2)





padj_threshold <-  0.05
lfc_threshold <- 1.5

# Add significance categories
all_res$sig_category <- "Not Significant"
all_res$sig_category[all_res$padj < padj_threshold & all_res$log2FoldChange > lfc_threshold] <- "Upregulated"
all_res$sig_category[all_res$padj < padj_threshold & all_res$log2FoldChange < -lfc_threshold] <- "Downregulated"

# Calculate gene counts
total_genes <- nrow(all_res)
up_genes <- sum(all_res$sig_category == "Upregulated", na.rm = TRUE)
down_genes <- sum(all_res$sig_category == "Downregulated", na.rm = TRUE)
sig_genes <- up_genes + down_genes

# Create gene count labels
stat_label <-sprintf("Total: %d\nUp: %d\nDown: %d\nPadj < 0.05 and log2FC > 1.5", 
                     total_genes, up_genes, down_genes)


sig_genes <- subset(all_res, padj < 0.05 & abs(log2FoldChange) > 1.5)
top_genes <- sig_genes[order(sig_genes$padj)[1:min(100, nrow(sig_genes))], ]
top_genes <- as.list(top_genes$X)

p2 <- ggplot(all_res, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(colour = "grey", alpha = 0.1) + # Colour all points grey
  geom_point(data = all_res[all_res$log2FoldChange < 0 & all_res$padj <= 0.05,], 
             aes(x = log2FoldChange, y = -log10(padj)), colour = "purple", alpha = 0.1) + # Colour downreg, significant in blue
  geom_point(data = all_res[all_res$log2FoldChange > 0 & all_res$padj <= 0.05,], 
             aes(x = log2FoldChange, y = -log10(padj)), colour = "orange", alpha = 0.1) + # Colour upreg, significant in red
  geom_text_repel(data = all_res[all_res$gene_symbol %in% top_genes & abs(all_res$log2FoldChange) > 0,], 
                  aes(x = log2FoldChange, y = -log10(padj), label = gene_symbol),
                  max.overlaps = Inf,  # Allow all labels to be displayed
                  box.padding = 0.1,     # Increase spacing between labels
                  point.padding = 0.1,  # Increase spacing from points
                  min.segment.length = 0.1,  # Minimum line segment length
                  force = 0.1,          # Increase force of repulsion
                  seed = 123,          # Set seed for reproducibility
                  size = 2,            # Smaller text size
                  segment.size = 0.1,  # Thinner connection lines
                  segment.alpha = 0.8) +  # Slightly transparent line
  labs(title = "ICM: Ulixertinib vs DMSO", subtitle = stat_label)+
  theme_classic() + theme(aspect.ratio = 0.75)

plot(p2)
ggsave("ICM_volcano_top100padj.pdf", p2)



#===============================================================================
### PrE
#===============================================================================

# Load PrE all_res
all_res <- read.csv("PrE_Ulix_vs_DMSO_allDEGs.csv")

#===============================================================================
# Volcano plot of all DEGs
#===============================================================================

#===============================================================================
# Volcano plot EPI DEGs with genes of interest
#===============================================================================

PrE_GOI <- read.csv("PrE_genelist.csv", header = FALSE)
genes <- as.list(PrE_GOI$V1)


p2 <- ggplot(all_res, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(colour = "grey", alpha = 0.1) + # Colour all points grey
  geom_point(data = all_res[all_res$log2FoldChange < 0 & all_res$padj <= 0.05,], 
             aes(x = log2FoldChange, y = -log10(padj)), colour = "purple", alpha = 0.1) + # Colour downreg, significant in blue
  geom_point(data = all_res[all_res$log2FoldChange > 0 & all_res$padj <= 0.05,], 
             aes(x = log2FoldChange, y = -log10(padj)), colour = "orange", alpha = 0.1) + # Colour upreg, significant in red
  geom_text_repel(data = all_res[all_res$gene_symbol %in% genes & all_res$log2FoldChange >0,], 
                  aes(x = log2FoldChange, y = -log10(padj), label = gene_symbol), 
                  min.segment.length = unit(0, 'lines'), nudge_y = 1, nudge_x = 2, size = 3, fontface = "italic", direction = "y") +
  geom_text_repel(data = all_res[all_res$gene_symbol %in% genes & all_res$log2FoldChange <0,], 
                  aes(x = log2FoldChange, y = -log10(padj), label = gene_symbol), 
                  min.segment.length = unit(0, 'lines'), nudge_y = 1, nudge_x = -2, size = 3, fontface = "italic", direction = "y") +
  labs(title = "Ulixertinib vs DMSO", subtitle = "DESeq2 padj < 0.05")+
  theme_classic() + theme(aspect.ratio = 1)

plot(p2)
ggsave("PrE_volcano_GOI.pdf", p2)





padj_threshold <-  0.05
lfc_threshold <- 1.5

# Add significance categories
all_res$sig_category <- "Not Significant"
all_res$sig_category[all_res$padj < padj_threshold & all_res$log2FoldChange > lfc_threshold] <- "Upregulated"
all_res$sig_category[all_res$padj < padj_threshold & all_res$log2FoldChange < -lfc_threshold] <- "Downregulated"

# Calculate gene counts
total_genes <- nrow(all_res)
up_genes <- sum(all_res$sig_category == "Upregulated", na.rm = TRUE)
down_genes <- sum(all_res$sig_category == "Downregulated", na.rm = TRUE)
sig_genes <- up_genes + down_genes

# Create gene count labels
stat_label <-sprintf("Total: %d\nUp: %d\nDown: %d\nPadj < 0.05 and log2FC > 1.5", 
                     total_genes, up_genes, down_genes)


sig_genes <- subset(all_res, padj < 0.05 & abs(log2FoldChange) > 1.5)
top_genes <- sig_genes[order(sig_genes$padj)[1:min(100, nrow(sig_genes))], ]
top_genes <- as.list(top_genes$X)

p2 <- ggplot(all_res, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(colour = "grey", alpha = 0.1) + # Colour all points grey
  geom_point(data = all_res[all_res$log2FoldChange < 0 & all_res$padj <= 0.05,], 
             aes(x = log2FoldChange, y = -log10(padj)), colour = "purple", alpha = 0.1) + # Colour downreg, significant in blue
  geom_point(data = all_res[all_res$log2FoldChange > 0 & all_res$padj <= 0.05,], 
             aes(x = log2FoldChange, y = -log10(padj)), colour = "orange", alpha = 0.1) + # Colour upreg, significant in red
  geom_text_repel(data = all_res[all_res$gene_symbol %in% top_genes & abs(all_res$log2FoldChange) > 0,], 
                  aes(x = log2FoldChange, y = -log10(padj), label = gene_symbol),
                  max.overlaps = Inf,  # Allow all labels to be displayed
                  box.padding = 0.1,     # Increase spacing between labels
                  point.padding = 0.1,  # Increase spacing from points
                  min.segment.length = 0.1,  # Minimum line segment length
                  force = 0.1,          # Increase force of repulsion
                  seed = 123,          # Set seed for reproducibility
                  size = 2,            # Smaller text size
                  segment.size = 0.1,  # Thinner connection lines
                  segment.alpha = 0.8) +  # Slightly transparent line
  labs(title = "PrE: Ulixertinib vs DMSO", subtitle = stat_label)+
  theme_classic() + theme(aspect.ratio = 0.75)

plot(p2)
ggsave("PrE_volcano_top100padj.pdf", p2)








#===============================================================================
# Minimalist Epiblast Volcano plot for Fig. 3D
#===============================================================================
all_res <- read.csv("EPI_Ulix_vs_DMSO_allDEGs.csv")

padj_threshold <-  0.05
lfc_threshold <- 1.5

sig_genes <- subset(all_res, padj < 0.05 & abs(log2FoldChange) > 1.5)
top_genes <- sig_genes[order(sig_genes$padj)[1:min(50, nrow(sig_genes))], ]
top_genes <- as.list(top_genes$X)


# Add significance categories
all_res$sig_category <- "NS"
all_res$sig_category[all_res$padj < padj_threshold & all_res$log2FoldChange > lfc_threshold] <- "Up"
all_res$sig_category[all_res$padj < padj_threshold & all_res$log2FoldChange < -lfc_threshold] <- "Down"

# Calculate gene counts
total_genes <- nrow(all_res)
up_genes <- sum(all_res$sig_category == "Up", na.rm = TRUE)
down_genes <- sum(all_res$sig_category == "Down", na.rm = TRUE)
sig_genes <- up_genes + down_genes

p2 <- ggplot(all_res, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = sig_category), alpha = 0.5, size = 3) + # Colour all points grey
  # geom_point(data = all_res[all_res$log2FoldChange < -1.5 & all_res$padj < 0.05,], 
  #            aes(x = log2FoldChange, y = -log10(padj)), colour = "purple", alpha = 0.1, size = 6) + # Colour downreg, significant in blue
  # geom_point(data = all_res[all_res$log2FoldChange > 1.5 & all_res$padj < 0.05,], 
  #            aes(x = log2FoldChange, y = -log10(padj)), colour = "orange", alpha = 0.1, size = 6) + # Colour upreg, significant in red
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black", alpha = 0.5) + 
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "black",  alpha = 0.5) + 
  # geom_text_repel(data = all_res[all_res$gene_symbol %in% top_genes & abs(all_res$log2FoldChange) > 0,], 
  #                 aes(x = log2FoldChange, y = -log10(padj), label = gene_symbol),
  #                 max.overlaps = Inf,  
  #                 box.padding = 0.1,    
  #                 point.padding = 0.1, 
  #                 min.segment.length = 0.1, 
  #                 force = 0.1,        
  #                 seed = 123,          
  #                 size = 2,           
  #                 segment.size = 0.1,  
  #                 segment.alpha = 0.8) +  
annotate("text", x = -10, y = max(-log10(all_res$padj), na.rm = TRUE) - 1, 
         label = paste0("Down: ", down_genes), hjust = 0, size = 8) +
  annotate("text", x = 10, y = max(-log10(all_res$padj), na.rm = TRUE) - 1, 
           label = paste0("Up: ", up_genes), hjust = 1, size = 8) +
#  labs(title = "Epiblast: ERKi vs Control") +
  theme_classic() + theme(text = element_text(size = 32), aspect.ratio = 0.5, 
                          legend.title = element_blank()) + 
  scale_color_manual(values = c("NS" = "grey", 
                                "Up" = "orange", 
                                "Down" = "purple")) 
plot(p2)
ggsave("volcano_epi.png", p2, width = 4000, height = 2000, dpi = 300, units = "px", scale = 1)


#===============================================================================
# Minimalist Hypoblast Volcano plot for SUpplemental Fig. 3D
#===============================================================================
all_res <- read.csv("PrE_Ulix_vs_DMSO_allDEGs.csv")

padj_threshold <-  0.05
lfc_threshold <- 1.5

sig_genes <- subset(all_res, padj < 0.05 & abs(log2FoldChange) > 1.5)
top_genes <- sig_genes[order(sig_genes$padj)[1:min(50, nrow(sig_genes))], ]
top_genes <- as.list(top_genes$X)


# Add significance categories
all_res$sig_category <- "NS"
all_res$sig_category[all_res$padj < padj_threshold & all_res$log2FoldChange > lfc_threshold] <- "Up"
all_res$sig_category[all_res$padj < padj_threshold & all_res$log2FoldChange < -lfc_threshold] <- "Down"

# Calculate gene counts
total_genes <- nrow(all_res)
up_genes <- sum(all_res$sig_category == "Up", na.rm = TRUE)
down_genes <- sum(all_res$sig_category == "Down", na.rm = TRUE)
sig_genes <- up_genes + down_genes

p2 <- ggplot(all_res, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = sig_category), alpha = 0.5, size = 3) + # Colour all points grey
  # geom_point(data = all_res[all_res$log2FoldChange < -1.5 & all_res$padj < 0.05,], 
  #            aes(x = log2FoldChange, y = -log10(padj)), colour = "purple", alpha = 0.1, size = 6) + # Colour downreg, significant in blue
  # geom_point(data = all_res[all_res$log2FoldChange > 1.5 & all_res$padj < 0.05,], 
  #            aes(x = log2FoldChange, y = -log10(padj)), colour = "orange", alpha = 0.1, size = 6) + # Colour upreg, significant in red
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black", alpha = 0.5) + 
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "black",  alpha = 0.5) + 
  # geom_text_repel(data = all_res[all_res$gene_symbol %in% top_genes & abs(all_res$log2FoldChange) > 0,], 
  #                 aes(x = log2FoldChange, y = -log10(padj), label = gene_symbol),
  #                 max.overlaps = Inf,  
  #                 box.padding = 0.1,    
  #                 point.padding = 0.1, 
  #                 min.segment.length = 0.1, 
  #                 force = 0.1,        
  #                 seed = 123,          
  #                 size = 2,           
  #                 segment.size = 0.1,  
  #                 segment.alpha = 0.8) +  
annotate("text", x = -10, y = max(-log10(all_res$padj), na.rm = TRUE) - 1, 
         label = paste0(down_genes), hjust = 0, size = 8) +
  annotate("text", x = 10, y = max(-log10(all_res$padj), na.rm = TRUE) - 1, 
           label = paste0(up_genes), hjust = 1, size = 8) +
  labs(title = "Hypoblast: ERKi vs Control") +
  theme_classic() + theme(text = element_text(size = 32), aspect.ratio = 1, 
                          legend.title = element_blank()) + 
  scale_color_manual(values = c("NS" = "grey", 
                                "Up" = "orange", 
                                "Down" = "purple")) 
plot(p2)
ggsave("volcano_PrE.png", p2, width = 3000, height = 3000, dpi = 300, units = "px", scale = 1)


#===============================================================================
# Minimalist TE Volcano plot for Supplemental Fig. 3D
#===============================================================================
all_res <- read.csv("TE_Ulix_vs_DMSO_allDEGs.csv")

padj_threshold <-  0.05
lfc_threshold <- 1.5

# Create gene count labels
stat_label <-sprintf("Total: %d\nUp: %d\nDown: %d\nPadj < 0.05 and log2FC > 1.5", 
                     total_genes, up_genes, down_genes)


sig_genes <- subset(all_res, padj < 0.05 & abs(log2FoldChange) > 1.5)
top_genes <- sig_genes[order(sig_genes$padj)[1:min(50, nrow(sig_genes))], ]
top_genes <- as.list(top_genes$X)


# Add significance categories
all_res$sig_category <- "NS"
all_res$sig_category[all_res$padj < padj_threshold & all_res$log2FoldChange > lfc_threshold] <- "Up"
all_res$sig_category[all_res$padj < padj_threshold & all_res$log2FoldChange < -lfc_threshold] <- "Down"


# Calculate gene counts
total_genes <- nrow(all_res)
up_genes <- sum(all_res$sig_category == "Up", na.rm = TRUE)
down_genes <- sum(all_res$sig_category == "Down", na.rm = TRUE)
sig_genes <- up_genes + down_genes

p2 <- ggplot(all_res, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = sig_category), alpha = 0.5, size = 3) + # Colour all points grey
  # geom_point(data = all_res[all_res$log2FoldChange < -1.5 & all_res$padj < 0.05,], 
  #            aes(x = log2FoldChange, y = -log10(padj)), colour = "purple", alpha = 0.1, size = 6) + # Colour downreg, significant in blue
  # geom_point(data = all_res[all_res$log2FoldChange > 1.5 & all_res$padj < 0.05,], 
  #            aes(x = log2FoldChange, y = -log10(padj)), colour = "orange", alpha = 0.1, size = 6) + # Colour upreg, significant in red
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black", alpha = 0.5) + 
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "black",  alpha = 0.5) + 
  # geom_text_repel(data = all_res[all_res$gene_symbol %in% top_genes & abs(all_res$log2FoldChange) > 0,], 
  #                 aes(x = log2FoldChange, y = -log10(padj), label = gene_symbol),
  #                 max.overlaps = Inf,  
  #                 box.padding = 0.1,    
  #                 point.padding = 0.1, 
  #                 min.segment.length = 0.1, 
  #                 force = 0.1,        
  #                 seed = 123,          
  #                 size = 2,           
  #                 segment.size = 0.1,  
  #                 segment.alpha = 0.8) +  
annotate("text", x = -10, y = max(-log10(all_res$padj), na.rm = TRUE) - 1, 
         label = paste0(down_genes), hjust = 0, size = 8) +
  annotate("text", x = 10, y = max(-log10(all_res$padj), na.rm = TRUE) - 1, 
           label = paste0(up_genes), hjust = 1, size = 8) +
  labs(title = "TE: ERKi vs Control") +
  theme_classic() + theme(text = element_text(size = 32), aspect.ratio = 1, 
                          legend.title = element_blank()) + 
  scale_color_manual(values = c("NS" = "grey", 
                                "Up" = "orange", 
                                "Down" = "purple")) 
plot(p2)
ggsave("volcano_TE.png", p2, width = 3000, height = 3000, dpi = 300, units = "px", scale = 1)


