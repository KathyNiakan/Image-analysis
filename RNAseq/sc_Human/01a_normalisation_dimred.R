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

# Set work directory
setwd("~/Documents/Collaborations/Brickman/61224/11_niakan_adata_combined")

# Gene filtered raw counts generated in 00_gene_summary.filter.R
load("sce_filtered.RData")

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

p1 <- ggplot(sce_pca, aes(PC1, PC2, color=treatment)) + 
  geom_point(size = 2) +  theme_classic()
print(p1)
p1 <- ggplot(sce_pca, aes(PC1, PC2, color=prediction)) + 
  geom_point(size = 2) +  theme_classic()
print(p1)

# Colour pallete
cols <- c("Inner Cell Mass" = "#7F7F7F",
          "Epiblast_6.0" = "#00FF00", 
          "Epiblast_7.0" = "#00BF00", 
          "Late epiblast" = "#007F00", 
          "Primitive Endoderm" = "#FF00FF",
          "Trophectoderm_6.0" = "#0000FF",
          "Trophectoderm_7.0" = "#0000BF",
          "Trophectoderm_8.0" = "#00007F",
          "Trophectoderm_9.0" = "#00003F")

p1 <- ggplot(sce_pca, aes(PC1, PC2, color=prediction, shape=treatment)) + 
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
    geom_point(size=4, alpha = 0.8) +
    scale_color_viridis_c() +
    theme_classic() +
    labs(title=gene_symbol, color="Log Expression")
}

# Create multi-panel plot for all genes
hvg_genes <- c("NANOG", "SOX2", "KLF17", "TDGF1", "PDGFRA", "GATA6", 
               "GATA4", "SOX17", "GATA2", "GATA3", "KRT18", "TEAD3")

# Create plot list
plot_list <- lapply(hvg_genes, function(gene) plot_gene_pca(sce, gene))

# Combine plots using patchwork
combined_plot <- wrap_plots(plot_list, ncol = 2)

# Print and save the combined plot
print(combined_plot)
ggsave("241219_gene_expression_pca.pdf", combined_plot, width = 15, height = 20)


# ==============================================================================
# Plot clustering and gene expression on PCA
# ==============================================================================

# Create plot list for PCA gene expression
plot_list_pca <- lapply(hvg_genes, function(gene) plot_gene_pca(sce, gene))

# Create the combined layout
layout <- (p1) / wrap_plots(plot_list_pca, ncol = 2)

# Add titles to the cluster and prediction plots
p1 <- p1 + ggtitle("Predictions")

# Combine all plots
final_plot <- layout +
  plot_layout(heights = c(1, 3)) &
  theme(plot.title = element_text(hjust = 0.5))
plot(final_plot)
# Save as PDF
ggsave("241219_combined_pca_analysis.pdf", final_plot, width = 16, height = 30)

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
# Plot UMAP with cluster annotation / prediction from Brickman lab
# ==============================================================================

p1 <- ggplot(sce_umap, aes(UMAP1, UMAP2, color=treatment)) + 
  geom_point(size = 2) +  theme_classic()
print(p1)

p1 <- ggplot(sce_umap, aes(UMAP1, UMAP2, color=prediction)) + 
  geom_point(size = 2) +  theme_classic()
print(p1)

p1 <- ggplot(sce_umap, aes(UMAP1, UMAP2, color=prediction, shape=treatment)) + 
  geom_point(size = 3, alpha = 0.6) +  theme_classic() + scale_color_manual(values = cols) 
print(p1)

# ==============================================================================
# Re-annotate clusters based on new UMAP clustering
# ==============================================================================

# Perform k-means clustering on UMAP coordinates (or PCA coordinates)
set.seed(1234)
num_clusters <- 9  # Specify the number of clusters

# Perform k-means clustering on the UMAP results
kmeans_res <- kmeans(sce_umap[, c("UMAP1", "UMAP2")], centers = num_clusters)
sce_umap$cluster <- factor(kmeans_res$cluster)

# ==============================================================================
# Plot UMAP with new cluster annotation
# ==============================================================================

p1 <- ggplot(sce_umap, aes(x = UMAP1, y = UMAP2, color = cluster, shape=treatment)) +
  geom_point(size = 2, alpha = 0.8) +  # Adjust size and transparency
  theme_classic() + 
  labs(title = "UMAP with Clusters", color = "Cluster")  # Add title and label
print(p1)

# ==============================================================================
# Manual re-annotation of k-means clusters
# ==============================================================================

sce_umap$cluster_label <- as.character(sce_umap$cluster)
sce_umap$cluster_label[sce_umap$cluster == 1] <- "TE"
sce_umap$cluster_label[sce_umap$cluster == 2] <- "EPI"
sce_umap$cluster_label[sce_umap$cluster == 3] <- "TE"
sce_umap$cluster_label[sce_umap$cluster == 4] <- "TE"
sce_umap$cluster_label[sce_umap$cluster == 5] <- "TE"
sce_umap$cluster_label[sce_umap$cluster == 6] <- "TE"
sce_umap$cluster_label[sce_umap$cluster == 7] <- "TE"
sce_umap$cluster_label[sce_umap$cluster == 8] <- "TE"
sce_umap$cluster_label[sce_umap$cluster == 9] <- "PrE"

# Note the multiple clusters are in the TE

# ==============================================================================
# Plot new and old clustering on UMAP, with Log gene expression
# ==============================================================================

# New clustering colour pallete
clust <- c("EPI" = "#00FF00", 
           "PrE" = "#FF00FF",
           "TE" = "#0000FF")

# Plot the UMAP with the manual annotations
p1 <- ggplot(sce_umap, aes(x = UMAP1, y = UMAP2, color = cluster_label, shape = treatment)) +
  geom_point(size = 4, alpha = 0.8) + scale_color_manual(values = clust) +
  theme_classic()
print(p1)

# Plot the UMAP with the manual annotations
p2 <- ggplot(sce_umap, aes(x = UMAP1, y = UMAP2, color = prediction, shape = treatment)) +
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
    treatment = sce$treatment
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
ggsave("241219_gene_expression_umap.pdf", combined_plot_umap, width = 15, height = 20)


# Create plot list for UMAP gene expression
plot_list_umap <- lapply(hvg_genes, function(gene) plot_gene_umap(sce, gene))

# Create the combined layout
layout <- (p1 + p2) / wrap_plots(plot_list_umap, ncol = 2)

# Add titles to the cluster and prediction plots
p1 <- p1 + ggtitle("Clusters")
p2 <- p2 + ggtitle("Predictions")

# Combine all plots
final_plot <- layout +
  plot_layout(heights = c(1, 3)) &
  theme(plot.title = element_text(hjust = 0.5))
print(final_plot)
# Save as PDF
ggsave("241219_combined_umap_analysis.pdf", final_plot, width = 16, height = 30)


# ==============================================================================
# UMAP, with log gene expression. Figure 3A, B and C, Supp Fig 3A
# ==============================================================================


treat <- c("DMSO" = 21, 
           "Ulix" = 24)

# Plot the UMAP with the manual annotations
p1 <- ggplot(sce_umap, aes(x = UMAP1, y = UMAP2, fill = cluster_label, shape = treatment)) +
  geom_point(size = 8, alpha = 0.2, stroke = 0.5) + scale_fill_manual(values = clust) + 
  scale_shape_manual(values = treat) +
theme_classic() + theme(aspect.ratio = 1,
                        text = element_text(size = 32))

print(p1)
ggsave("UMAP_lineage.png", p1, width = 2350, height = 2350, dpi = 300, units = "px", scale = 1)
ggsave("UMAP_lineage.pdf", p1, width = 2350, height = 2350, dpi = 300, units = "px", scale = 1)


coltreat <- c("DMSO" = "black", 
           "Ulix" = "orange")
clustshape <- c("EPI" = 21, 
           "PrE" = 21,
           "TE" = 21)

# Plot the UMAP with treatment
p1 <- ggplot(sce_umap, aes(x = UMAP1, y = UMAP2, fill = treatment, shape = cluster_label)) +
  geom_point(size = 8, alpha = 0.2, stroke = 0.5) + scale_fill_manual(values = coltreat) + scale_shape_manual(values = clustshape) +
  theme_classic() + theme(aspect.ratio = 1,
                          text = element_text(size = 32))

print(p1)
ggsave("UMAP_treat.png", p1, width = 2350, height = 2350, dpi = 300, units = "px", scale = 1)


main_fig_genes <- c("NANOG", "SOX2", 
               "GATA4", "SOX17")

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
    treatment = sce$treatment
  )
  
  plot_df <- plot_df[order(plot_df$expression), ]
  
  ggplot(plot_df, aes(UMAP1, UMAP2, color=expression)) + 
    geom_point(size=8, alpha = 0.8) +
    scale_color_viridis_c() +
    theme_classic() + theme(aspect.ratio = 1,
                            text = element_text(size = 32),
                            legend.title = element_text(face = "italic")) + 
    labs(color = gene_symbol)
}

# Create plot list for UMAP
plot_list_umap <- lapply(main_fig_genes, function(gene) plot_gene_umap(sce, gene))

# Combine UMAP plots
combined_plot_umap <- wrap_plots(plot_list_umap, ncol = 2)

# Print and save
print(combined_plot_umap)
ggsave("UMAP_genes.png", combined_plot_umap, width = 6000, height = 3000, dpi = 300, units = "px", scale = 1)



# Create multi-panel plot for all genes
supp_genes <- c( "KLF17", "TDGF1", "PDGFRA", "GATA6", 
               "GATA2", "GATA3", "KRT18", "TEAD3")

# Create plot list for UMAP
suppplot_list_umap <- lapply(supp_genes, function(gene) plot_gene_umap(sce, gene))

# Combine UMAP plots
supp_plot_umap <- wrap_plots(suppplot_list_umap, ncol = 4)

# Print and save
print(combined_plot_umap)
ggsave("UMAP_suppgenes.png", supp_plot_umap, width = 12000, height = 3000, dpi = 300, units = "px", scale = 1)


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
    treatment = sce$treatment
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
  group_by(treatment, cluster_label) %>%
  summarize(N = n()) %>%
  ungroup() %>%
  complete(treatment, cluster_label,
           fill = list(N = 0))

g <- ICM_count %>%
  mutate(cluster_label = factor(cluster_label, levels=c("EPI", "PrE"))) %>%
  ggplot(aes(x = treatment, y = N, fill = cluster_label))
g <- g + geom_bar(position="fill", stat="identity")
g <- g + scale_fill_manual(values = clust) + theme_classic() + 
  scale_y_continuous(labels = scales::percent) +
  ylab("% of ICM Cells") + theme(text = element_text(size = 8)) + theme(legend.position="none",
                                                                        axis.title.x = element_blank(),
                                                                        axis.text.x = element_text(angle = 45, hjust=1),
                                                                        aspect.ratio = 2)
print(g)
ggsave("ICMlineage_stacked.png", g, width = 350, height = 500, dpi = 300, units = "px", scale = 1)



# ==============================================================================
# Counts similar to stacked ICM for IF. Brickman lab annotation.
# ==============================================================================


# Colour pallete
matched <- c("Inner Cell Mass" = "#F6C445",
          "Epiblast_6.0" = "#D6B2CA", 
          "Epiblast_7.0" = "#C38DB1", 
          "Late epiblast" = "#AA5C8F", 
          "Primitive Endoderm" = "#D05B61",
          "Trophectoderm_6.0" = "#BDD4EB",
          "Trophectoderm_7.0" = "#ACC9E6",
          "Trophectoderm_8.0" = "#9CBFE2",
          "Trophectoderm_9.0" = "#8BB4DD")



ICM_count_B <- cell_metadata_df %>% 
  group_by(treatment, prediction) %>%
  summarize(N = n()) %>%
  ungroup() %>%
  complete(treatment, prediction,
           fill = list(N = 0))

g <- ICM_count_B %>%
  mutate(prediction = factor(prediction, levels=c("Inner Cell Mass",
          "Epiblast_6.0", 
          "Epiblast_7.0", 
          "Late epiblast", 
          "Primitive Endoderm",
          "Trophectoderm_6.0",
          "Trophectoderm_7.0",
          "Trophectoderm_8.0",
          "Trophectoderm_9.0") )) %>%
  ggplot(aes(x = treatment, y = N, fill = prediction))
g <- g + geom_bar(position="fill", stat="identity")
g <- g + scale_fill_manual(values = matched) + theme_classic(base_size = 22) + 
  scale_y_continuous(labels = scales::percent) +
  ylab("% of Cells") + theme(text = element_text(size = 18)) + theme(legend.position = "none",
                                                                     axis.title.x = element_blank(),
                                                                        axis.text.x = element_text(angle = 45, hjust=1),
                                                                        aspect.ratio = 2)
print(g)
ggsave("lineage_stacked_prediction_matched.png", g, width = 1200, height = 1200, dpi = 300, units = "px", scale = 1)


g <- ICM_count_B %>%   filter(!grepl("Trophectoderm", prediction)) %>% 
  mutate(prediction = factor(prediction, levels=c("Inner Cell Mass",
                                                  "Epiblast_6.0", 
                                                  "Epiblast_7.0", 
                                                  "Late epiblast", 
                                                  "Primitive Endoderm") )) %>%
  ggplot(aes(x = treatment, y = N, fill = prediction))
g <- g + geom_bar(position="fill", stat="identity", alpha = 0.6, color = "black")
g <- g + scale_fill_manual(values = cols) + theme_classic() + 
  scale_y_continuous(labels = scales::percent) +
  ylab("% of Cells") + theme(text = element_text(size = 8)) + theme(axis.title.x = element_blank(),
                                                                    axis.text.x = element_text(angle = 45, hjust=1),
                                                                    aspect.ratio = 2)
print(g)
ggsave("ICMlineage_stacked_prediction.png", g, width = 350, height = 500, dpi = 300, units = "px", scale = 1)


g <- ICM_count_B %>%   filter(grepl("piblast", prediction)) %>% 
  mutate(prediction = factor(prediction, levels=c("Epiblast_6.0", 
                                                  "Epiblast_7.0", 
                                                  "Late epiblast") )) %>%
  ggplot(aes(x = treatment, y = N, fill = prediction))
g <- g + geom_bar(position="fill", stat="identity", alpha = 0.6, color = "black")
g <- g + scale_fill_manual(values = cols) + theme_classic() + 
  scale_y_continuous(labels = scales::percent) +
  ylab("% of Cells") + theme(text = element_text(size = 8)) + theme(axis.title.x = element_blank(),
                                                                    axis.text.x = element_text(angle = 45, hjust=1),
                                                                    aspect.ratio = 2)
print(g)
ggsave("EPIlineage_stacked_prediction.png", g, width = 350, height = 500, dpi = 300, units = "px", scale = 1)


g <- ICM_count_B %>%   filter(grepl("piblast", prediction)) %>% 
  mutate(prediction = factor(prediction, levels=c("Epiblast_6.0", 
                                                  "Epiblast_7.0", 
                                                  "Late epiblast") )) %>%
  ggplot(aes(x = treatment, y = N, fill = prediction))
g <- g + geom_bar(stat="identity", alpha = 0.6, color = "black")
g <- g + scale_fill_manual(values = cols) + theme_classic() + 
  scale_y_continuous() +
  ylab("Number of Cells") + theme(text = element_text(size = 8)) + theme(axis.title.x = element_blank(),
                                                                    axis.text.x = element_text(angle = 45, hjust=1),
                                                                    aspect.ratio = 2)
print(g)



g <- ICM_count_B %>% filter(!grepl("Trophectoderm", prediction)) %>% 
  mutate(prediction = factor(prediction, levels=c("Inner Cell Mass",
                                                  "Epiblast_6.0", 
                                                  "Epiblast_7.0", 
                                                  "Late epiblast", 
                                                  "Primitive Endoderm") )) %>%
  ggplot(aes(x = treatment, y = N, fill = prediction))
g <- g + geom_bar(stat="identity", alpha = 0.6, color = "black")
g <- g + scale_fill_manual(values = cols) + theme_classic() + 
  scale_y_continuous() +
  ylab("Number of ICM Cells") + theme(text = element_text(size = 8)) + theme(axis.title.x = element_blank(),
                                                                         axis.text.x = element_text(angle = 45, hjust=1),
                                                                         aspect.ratio = 2)
print(g)