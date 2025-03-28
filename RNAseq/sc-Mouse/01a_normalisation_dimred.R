# ==============================================================================
# Load packages, set wd
# ==============================================================================

library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)
library(patchwork)
library(dplyr)
library(annotables)
library(biomaRt)
library(SummarizedExperiment)
library(tidyr)

# Set work directory
setwd("~/Documents/Collaborations/Brickman/121224")

# Gene filtered raw counts generated in 00_gene_summary.filter.R
load("sce_filtered.RData")


# ==============================================================================
# Sample filtering
# ==============================================================================

# After checking gene expression there is clearly one TE outlier among these genes
# Residual after immunosurgery

TE_cell <- c("SIM5111A80_SIM5111A80")
sce <- sce[, !(colnames(sce) %in% TE_cell)]


# ==============================================================================
# Normalisation and parameter set up for PCA
# ==============================================================================

# Size-factor normalisation
clust_sce <- quickCluster(sce, min.size = 10, d = 10, method = "igraph")
sce <- computeSumFactors(sce, cluster = clust_sce)
sce <- logNormCounts(sce)

# Parameters 
prop <- 0.1 # Proportion of highly variable genes (HVG) selected 10%
dims <- 25 # Number of components / dimensions

# Identify HVGs
dec <- modelGeneVar(sce)  # Model gene variance
hvg <- getTopHVGs(dec, prop = prop)  # Select top HVGs
rowData(sce)$is_hvg <- rownames(sce) %in% hvg

# ==============================================================================
# Run PCA
# ==============================================================================

# Run PCA on the HVGs
set.seed(1234)
sce <- runPCA(sce, ncomponents = dims, subset_row = hvg)

sce_pca <- reducedDim(sce, "PCA")
sce_pca <- sce_pca[,c("PC1", "PC2")]
meta <- data.frame(colData(sce))
sce_pca <- merge(sce_pca, meta, by = "row.names")
rownames(sce_pca) <- as.list(sce_pca$Row.names)
sce_pca$Row.names <- NULL

# ==============================================================================
# Plot PCAs
# ==============================================================================

p1 <- ggplot(sce_pca, aes(PC1, PC2, color=Treatment)) + 
  geom_point(size = 2) +  theme_classic()
print(p1)
p1 <- ggplot(sce_pca, aes(PC1, PC2, color=predictions)) + 
  geom_point(size = 2) +  theme_classic()
print(p1)


# Colour pallete
cols <- c("E3.5-ICM" = "#7F7F7F",
          "E4.5-EPI" = "#00FF00", 
          "E4.5-PrE" = "#FF00FF",
          "E3.5-TE" = "#0000FF")

p1 <- ggplot(sce_pca, aes(PC1, PC2, color=predictions, shape=Treatment)) + 
  geom_point(size = 5, alpha = 0.8) +  theme_classic() + scale_color_manual(values = cols) 
print(p1)



# ==============================================================================
# Plot gene expression on PCA
# ==============================================================================

plot_gene_pca <- function(sce, gene_symbol) {
  # Get ENSEMBL ID for the gene symbol
  ensembl_id <- rownames(sce)[which(rowData(sce)$gene_symbol == gene_symbol)]
  
  if(length(ensembl_id) == 0) {
    warning(paste("Gene", gene_symbol, "not found"))
    return(NULL)
  }
  
  # Get normalized expression
  expr <- logcounts(sce)[ensembl_id,]
  
  # Combine PCA coordinates with expression
  sce_pca <- data.frame(reducedDim(sce, "PCA")[,1:2])
  sce_pca$expression <- expr
  
  # Create plot
  ggplot(sce_pca, aes(PC1, PC2, color=expression)) + 
    geom_point(size=3, alpha = 0.8) +
    scale_color_viridis_c() +
    theme_classic() +
    labs(title=gene_symbol, color="Log Expression")
}

# Create multi-panel plot for all genes
hvg_genes <- c("Nanog", "Sox2", "Klf4", "Fgf4", 
               "Otx2", 
               "Pdgfra", "Gata6", "Fgfr2",
               "Gata4",  "Sox17", "Sox7"
               #, 
               #"Cdx2", "Gata3"
               )

# Create plot list
plot_list <- lapply(hvg_genes, function(gene) plot_gene_pca(sce, gene))

# Combine plots using patchwork
combined_plot <- wrap_plots(plot_list, ncol = 2)

# Print and save the combined plot
print(combined_plot)
ggsave("gene_expression_pca.pdf", combined_plot, width = 15, height = 20)


# ==============================================================================
# Plot clustering and gene expression on PCA
# ==============================================================================

# Create plot list for PCA gene expression
plot_list_pca <- lapply(hvg_genes, function(gene) plot_gene_pca(sce, gene))

# Create the combined layout
layout <- (p1) / wrap_plots(plot_list_pca, ncol = 2)

# Add titles to the cluster and predictions plots
p1 <- p1 + ggtitle("predictionss")

# Combine all plots
final_plot <- layout +
  plot_layout(heights = c(1, 3)) &
  theme(plot.title = element_text(hjust = 0.5))
plot(final_plot)
# Save as PDF
ggsave("combined_pca_analysis.pdf", final_plot, width = 16, height = 30)

# ==============================================================================
# UMAP
# ==============================================================================

# Run UMAP based on PCA analysis
sce <- runUMAP(sce, dimred = "PCA")
sce_umap <- reducedDim(sce, "UMAP")

sce_umap <- data.frame(sce_umap)
colnames(sce_umap) <- c("UMAP1", "UMAP2")

# Add sample info to UMAP data
meta <- data.frame(colData(sce))
sce_umap <- merge(sce_umap, meta, by = "row.names")
rownames(sce_umap) <- as.list(sce_umap$Row.names)
sce_umap$Row.names <- NULL

# ==============================================================================
# Plot UMAP with cluster annotation / predictions from Brickman lab
# ==============================================================================

p1 <- ggplot(sce_umap, aes(UMAP1, UMAP2, color=Treatment)) + 
  geom_point(size = 2) +  theme_classic()
print(p1)

p1 <- ggplot(sce_umap, aes(UMAP1, UMAP2, color=predictions)) + 
  geom_point(size = 2) +  theme_classic()
print(p1)

p1 <- ggplot(sce_umap, aes(UMAP1, UMAP2, color=predictions, shape=Treatment)) + 
  geom_point(size = 3, alpha = 0.6) +  theme_classic() + scale_color_manual(values = cols) 
print(p1)

# ==============================================================================
# Re-annotate clusters based on new UMAP clustering
# ==============================================================================

# Perform k-means clustering on UMAP coordinates (or PCA coordinates)
set.seed(1234)
num_clusters <- 2  # Specify the number of clusters

# Perform k-means clustering on the UMAP results
kmeans_res <- kmeans(sce_umap[, c("UMAP1", "UMAP2")], centers = num_clusters)
sce_umap$cluster <- factor(kmeans_res$cluster)

# ==============================================================================
# Plot UMAP with new cluster annotation
# ==============================================================================

p1 <- ggplot(sce_umap, aes(x = UMAP1, y = UMAP2, color = cluster, shape=Treatment)) +
  geom_point(size = 2, alpha = 0.8) +  # Adjust size and transparency
  theme_classic() + 
  labs(title = "UMAP with Clusters", color = "Cluster")  # Add title and label
print(p1)

# ==============================================================================
# Manual re-annotation of k-means clusters
# ==============================================================================

sce_umap$cluster_label <- as.character(sce_umap$cluster)
sce_umap$cluster_label[sce_umap$cluster == 1] <- "EPI"
sce_umap$cluster_label[sce_umap$cluster == 2] <- "PrE"


# ==============================================================================
# Plot new and old clustering on UMAP, with Log gene expression
# ==============================================================================

# New clustering colour pallete
clust <- c("EPI" = "#00FF00", 
           "PrE" = "#FF00FF"
           #,
          # "TE" = "#0000FF"
           )

# Plot the UMAP with the manual annotations
p1 <- ggplot(sce_umap, aes(x = UMAP1, y = UMAP2, color = cluster_label, shape = Treatment)) +
  geom_point(size = 4, alpha = 0.8) + scale_color_manual(values = clust) +
  theme_classic()
print(p1)

# Plot the UMAP with the manual annotations
p2 <- ggplot(sce_umap, aes(x = UMAP1, y = UMAP2, color = predictions, shape = Treatment)) +
  geom_point(size = 4, alpha = 0.8) + scale_color_manual(values = cols) +
  theme_classic()
print(p2)

# Plot gene expression on UMAP logCounts
plot_gene_umap <- function(sce, gene_symbol) {
  ensembl_id <- rownames(sce)[which(rowData(sce)$gene_symbol == gene_symbol)]
  
  if(length(ensembl_id) == 0) {
    warning(paste("Gene", gene_symbol, "not found"))
    return(NULL)
  }
  
  expr <- logcounts(sce)[ensembl_id,]
  umap_coords <- reducedDim(sce, "UMAP")
  
  plot_df <- data.frame(
    UMAP1 = umap_coords[,1],
    UMAP2 = umap_coords[,2],
    expression = expr,
    Treatment = sce$Treatment
  )
  
  ggplot(plot_df, aes(UMAP1, UMAP2, color=expression)) + 
    geom_point(size=4, alpha = 0.8) +
    scale_color_viridis_c() +
    theme_classic() +
    labs(title=gene_symbol, color="Log Expression")
}

# Create plot list for UMAP
plot_list_umap <- lapply(hvg_genes, function(gene) plot_gene_umap(sce, gene))

# Combine UMAP plots
combined_plot_umap <- wrap_plots(plot_list_umap, ncol = 2)

# Print and save
print(combined_plot_umap)
ggsave("gene_expression_umap.pdf", combined_plot_umap, width = 15, height = 20)

# Create plot list for UMAP gene expression
plot_list_umap <- lapply(hvg_genes, function(gene) plot_gene_umap(sce, gene))

# Create the combined layout
layout <- (p1 + p2) / wrap_plots(plot_list_umap, ncol = 2)

# Add titles to the cluster and predictions plots
p1 <- p1 + ggtitle("Clusters")
p2 <- p2 + ggtitle("predictionss")

# Combine all plots
final_plot <- layout +
  plot_layout(heights = c(1, 3)) &
  theme(plot.title = element_text(hjust = 0.5))
print(final_plot)
# Save as PDF
ggsave("combined_umap_analysis.pdf", final_plot, width = 16, height = 30)


# ==============================================================================
# UMAP, with log gene expression. Supp Fig 4
# ==============================================================================


treat <- c("DMSO" = 21, 
           "Ulixertinib" = 24)

# Plot the UMAP with the manual annotations
p1 <- ggplot(sce_umap, aes(x = UMAP1, y = UMAP2, fill = cluster_label, shape = Treatment)) +
  geom_point(size = 6, alpha = 0.6, stroke = 0.5) + scale_fill_manual(values = clust) + 
  scale_shape_manual(values = treat) +
theme_classic() + theme(aspect.ratio = 1,
                        text = element_text(size = 32))

print(p1)
ggsave("UMAP_lineage.png", p1, width = 2200, height = 2200, dpi = 300, units = "px", scale = 1)
ggsave("UMAP_lineage.pdf", p1, width = 2200, height = 2200, dpi = 300, units = "px", scale = 1)


coltreat <- c("DMSO" = "black", 
           "Ulixertinib" = "orange")
clustshape <- c("EPI" = 21, 
           "PrE" = 21,
           "TE" = 21)

# Plot the UMAP with Treatment
p1 <- ggplot(sce_umap, aes(x = UMAP1, y = UMAP2, fill = Treatment, shape = cluster_label)) +
  geom_point(size = 6, alpha = 0.6, stroke = 0.5) + scale_fill_manual(values = coltreat) + scale_shape_manual(values = clustshape) +
  theme_classic() + theme(aspect.ratio = 1,
                          text = element_text(size = 32))

print(p1)
ggsave("UMAP_treat.png", p1, width = 2200, height = 2200, dpi = 300, units = "px", scale = 1)


main_fig_genes <- c("Nanog", "Sox2", "Fgf4", 
               "Gata6", "Gata4", "Sox7")

# Plot gene expression on UMAP logCounts
plot_gene_umap <- function(sce, gene_symbol) {
  ensembl_id <- rownames(sce)[which(rowData(sce)$gene_symbol == gene_symbol)]
  
  if(length(ensembl_id) == 0) {
    warning(paste("Gene", gene_symbol, "not found"))
    return(NULL)
  }
  
  expr <- logcounts(sce)[ensembl_id,]
  umap_coords <- reducedDim(sce, "UMAP")
  
  plot_df <- data.frame(
    UMAP1 = umap_coords[,1],
    UMAP2 = umap_coords[,2],
    expression = expr,
    Treatment = sce$Treatment
  )
  
  plot_df <- plot_df[order(plot_df$expression), ]
  
  ggplot(plot_df, aes(UMAP1, UMAP2, color=expression)) + 
    geom_point(size=6, alpha = 0.8) +
    scale_color_viridis_c() +
    theme_classic() + theme(aspect.ratio = 1,
                            text = element_text(size = 32),
                            legend.title = element_text(face = "italic")) + 
    labs(color = gene_symbol)
}

# Create plot list for UMAP
plot_list_umap <- lapply(main_fig_genes, function(gene) plot_gene_umap(sce, gene))

# Combine UMAP plots
combined_plot_umap <- wrap_plots(plot_list_umap, ncol = 3)

# Print and save
print(combined_plot_umap)
ggsave("UMAP_genes.png", combined_plot_umap, width = 6000, height = 2800, dpi = 300, units = "px", scale = 1)


# ==============================================================================
# Plot new and old clustering on UMAP, with log TPM gene expression
# ==============================================================================

# Plot gene expression on UMAP - TPM
plot_gene_umap <- function(sce, gene_symbol) {
  ensembl_id <- rownames(sce)[which(rowData(sce)$gene_symbol == gene_symbol)]
  
  if(length(ensembl_id) == 0) {
    warning(paste("Gene", gene_symbol, "not found"))
    return(NULL)
  }
  
  expr <- assay(sce, "logtpm")[ensembl_id,]
  umap_coords <- reducedDim(sce, "UMAP")
  
  plot_df <- data.frame(
    UMAP1 = umap_coords[,1],
    UMAP2 = umap_coords[,2],
    expression = expr,
    Treatment = sce$Treatment
  )
  
  ggplot(plot_df, aes(UMAP1, UMAP2, color=expression)) + 
    geom_point(size=4, alpha = 0.8) +
    scale_color_viridis_c() +
    theme_classic() +
    labs(title=gene_symbol, color="Log TPM")
}

# Create plot list for UMAP
plot_list_umap <- lapply(hvg_genes, function(gene) plot_gene_umap(sce, gene))

# Combine UMAP plots
combined_plot_umap <- wrap_plots(plot_list_umap, ncol = 2)

# Print and save
print(combined_plot_umap)
ggsave("241219_gene_expression_umap_tpm.pdf", combined_plot_umap, width = 15, height = 20)

sce_umap <- sce_umap[colnames(sce), ]

# Add the cluster labels to the colData of the sce object
colData(sce)$cluster_label <- sce_umap$cluster_label

save(sce, file = "sce_filtered_cscluster.RData")

# Get metadata
cell_annotations <- colData(sce)

# Convert to a data frame
cell_metadata_df <- as.data.frame(cell_annotations)

# Save cell metadata as a CSV
write.csv(cell_metadata_df, "cell_metadata_cs.csv", row.names = TRUE)

# ==============================================================================
# Counts similar to stacked ICM for IF. Supp Fig. S4B
# ==============================================================================

ICM_count <- cell_metadata_df %>% filter(cluster_label != "TE") %>%
  group_by(Treatment, cluster_label) %>%
  summarize(N = n()) %>%
  ungroup() %>%
  complete(Treatment, cluster_label,
           fill = list(N = 0))

g <- ICM_count %>%
  mutate(cluster_label = factor(cluster_label, levels=c("EPI", "PrE"))) %>%
  ggplot(aes(x = Treatment, y = N, fill = cluster_label))
g <- g + geom_bar(position="fill", stat="identity")
g <- g + scale_fill_manual(values = clust) + theme_classic() + 
  scale_y_continuous(labels = scales::percent) +
  ylab("% of ICM Cells") + theme(text = element_text(size = 8)) + theme(legend.position="none",
                                                                        axis.title.x = element_blank(),
                                                                        axis.text.x = element_text(angle = 45, hjust=1),
                                                                        aspect.ratio = 2)
print(g)
ggsave("ICMlineage_stacked.png", g, width = 350, height = 500, dpi = 300, units = "px", scale = 1)





# Colour pallete
matched <- c("E3.5-ICM" = "#F8D06A",
             "E4.5-EPI" = "#B46F9C", 
             "E4.5-PrE" = "#D05B61")



ICM_count_B <- cell_metadata_df %>% 
  group_by(Treatment, predictions) %>%
  summarize(N = n()) %>%
  ungroup() %>%
  complete(Treatment, predictions,
           fill = list(N = 0))

g <- ICM_count_B %>%
  mutate(predictions = factor(predictions, levels=c("E3.5-ICM",
                                                    "E4.5-EPI", 
                                                    "E4.5-PrE") )) %>%
  ggplot(aes(x = Treatment, y = N, fill = predictions))
g <- g + geom_bar(position="fill", stat="identity")
g <- g + scale_fill_manual(values = matched) + theme_classic() + 
  scale_y_continuous(labels = scales::percent) +
  ylab("% of Cells") + theme(text = element_text(size = 8)) + theme(legend.position = "none",
                                                                     axis.title.x = element_blank(),
                                                                     axis.text.x = element_text(angle = 45, hjust=1),
                                                                     aspect.ratio = 2)
print(g)
ggsave("lineage_stacked_prediction_matched.png", g, width = 350, height = 500, dpi = 300, units = "px", scale = 1)


g <- ICM_count_B %>%   
  mutate(predictions = factor(predictions, levels=c("E3.5-ICM",
                                                    "E4.5-EPI", 
                                                    "E4.5-PrE") )) %>%
  ggplot(aes(x = Treatment, y = N, fill = predictions))
g <- g + geom_bar(stat="identity", alpha = 0.6, color = "black")
g <- g + scale_fill_manual(values = cols) + theme_classic() + 
  scale_y_continuous() +
  ylab("Number of Cells") + theme(text = element_text(size = 8)) + theme(axis.title.x = element_blank(),
                                                                         axis.text.x = element_text(angle = 45, hjust=1),
                                                                         aspect.ratio = 2)
print(g)



# Plot the UMAP with the manual annotations
p1 <- ggplot(sce_umap, aes(x = UMAP1, y = UMAP2, fill = predictions, shape = Treatment)) +
  geom_point(size = 6, alpha = 0.9, stroke = 0.5) + scale_fill_manual(values = matched) + 
  scale_shape_manual(values = treat) +
  theme_classic() + theme(aspect.ratio = 1,
                          text = element_text(size = 32))

print(p1)
ggsave("UMAP_predictions.png", p1, width = 2200, height = 2200, dpi = 300, units = "px", scale = 1)
ggsave("UMAP_predictions.pdf", p1, width = 2200, height = 2200, dpi = 300, units = "px", scale = 1)



sub <- c("sample", "cluster_label")
new_label <- cell_metadata_df[sub]
sce_pca <- left_join(sce_pca, new_label, by="sample")

p1 <- ggplot(sce_pca, aes(PC1, PC2, color=cluster_label, shape=Treatment)) + 
  geom_point(size = 5, alpha = 0.8) +  theme_classic() + scale_color_manual(values = clust) 
print(p1)


# Calculate the percentage of variance explained by each PC
pc_variance <- getExplanatoryPCs(sce, dimred = "PCA") [1:2, ]

# Plot the variance explained by PC1 and PC2
plotReducedDim(sce, dimred = "PCA", colour_by = "Treatment")

# The gene loadings are typically stored in `pca_loadings$rotation` (the eigenvectors)
pca_loadings <- reducedDims(sce)$PCA
loadings <- pca_loadings$rotation

# You can extract the top loading genes for PC1 and PC2
top_genes <- order(abs(loadings[, 1]), decreasing = TRUE)[1:10]  # top 10 genes for PC1
top_genes_pc2 <- order(abs(loadings[, 2]), decreasing = TRUE)[1:10]  # top 10 genes for PC2






recoded_meta <- cell_metadata_df 
recoded_meta$predictions <- recode_factor(recoded_meta$predictions, "E4.5-EPI" = "EPI")
recoded_meta$predictions <- recode_factor(recoded_meta$predictions, "E4.5-PrE" = "PrE")

counts <- recoded_meta %>% group_by(cluster_label, predictions) %>% summarize(N = n())

sum(counts[c(1,3),]$N) / sum(counts$N)