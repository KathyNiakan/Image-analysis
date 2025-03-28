# ==============================================================================
# Load packages, set wd and load data
# ==============================================================================
library(DESeq2)
library(readr)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tibble)
library(tidyr)
library(RColorBrewer)
library(EnhancedVolcano)
library(scran)
library(biomaRt)

# Set work directory
setwd("~/Documents/Collaborations/Brickman/61224/11_niakan_adata_combined")

# Gene filtered raw counts
load("sce_filtered_cscluster.RData")

# ==============================================================================
# Additional filtering steps
# ==============================================================================

# mitochondrial, ribosomal and duplicate genes already filtered out
# Now, filter out pseudo genes

# Get gene info
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_ids <- rownames(sce)

gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
filters = "ensembl_gene_id",
values = ensembl_ids,
mart = mart) %>%
  group_by(ensembl_gene_id) %>%
  summarize(
    gene_name = dplyr::first(external_gene_name),  
    gene_biotype = dplyr::first(gene_biotype)
  ) %>%
  as.data.frame()

gene_info$gene_name[gene_info$gene_name == ""] <- NA
no_name <- gene_info$ensembl_gene_id[is.na(gene_info$gene_name)]

# Filter out pseudogenes
pseudogenes <- gene_info %>% filter(grepl("pseudogene", gene_biotype))
sce <- sce[!(rownames(sce) %in% pseudogenes$ensembl_gene_id), ]
sce <- sce[!(rownames(sce) %in% no_name),]


# ==============================================================================
# Set up DEseq2 dataframe
# ==============================================================================
# Get gene annotations
gene_annotations <- as.data.frame(rowData(sce))

# Extract the counts matrix - keeping ensembl IDs as rownames
counts <- assay(sce, "counts")
counts_df <- as.data.frame(as.matrix(counts))
# No need to change rownames as they're already ensembl IDs
counts <- as.matrix(counts_df)

# Set up metadata
metadata <- as.data.frame(colData(sce))
# Modify the metadata to add ICM label
metadata$combined_cluster <- metadata$cluster_label
metadata$combined_cluster[metadata$cluster_label %in% c("EPI", "PrE")] <- "ICM"

# ==============================================================================
# Run DESeq2 analysis
# ==============================================================================

# Functions to create and run DESeq2 for a specific cluster, and pull out DEGs
source("functions/run_deseq.R")

# Run analysis for each cluster including ICM
clusters <- c("EPI", "TE", "PrE", "ICM")

# Get DEG analysis for each cluster
degs_list <- lapply(clusters, function(cluster) {
  message("Processing cluster: ", cluster)
  run_deseq_cluster(counts, metadata, cluster)
})
names(degs_list) <- clusters

# Access results for each cluster
epi_degs <- degs_list[["EPI"]]
te_degs <- degs_list[["TE"]]
pre_degs <- degs_list[["PrE"]]
icm_degs <- degs_list[["ICM"]]

# Function to add gene symbols to DEG results
add_gene_symbols <- function(deg_df) {
  # Add gene symbols column based on ensembl IDs
  deg_df$gene_symbol <- gene_annotations$gene_symbol[match(rownames(deg_df), rownames(gene_annotations))]
  return(deg_df)
}

# Get significant DEGs and add gene symbols
sig_degs <- lapply(degs_list, function(x) {
  add_gene_symbols(get_sig_degs(x))
})

# Export results with both ensembl IDs and gene symbols
write.csv(sig_degs[["EPI"]], "EPI_Ulix_vs_DMSO_DEGs.csv")
write.csv(sig_degs[["TE"]], "TE_Ulix_vs_DMSO_DEGs.csv")
write.csv(sig_degs[["PrE"]], "PrE_Ulix_vs_DMSO_DEGs.csv")
write.csv(sig_degs[["ICM"]], "ICM_Ulix_vs_DMSO_DEGs.csv")

# Get all DEGs and add gene symbols
all_degs <- lapply(degs_list, function(x) {
  add_gene_symbols(get_all_degs(x))
})

# Export all results
write.csv(all_degs[["EPI"]], "EPI_Ulix_vs_DMSO_allDEGs.csv")
write.csv(all_degs[["TE"]], "TE_Ulix_vs_DMSO_allDEGs.csv")
write.csv(all_degs[["PrE"]], "PrE_Ulix_vs_DMSO_allDEGs.csv")
write.csv(all_degs[["ICM"]], "ICM_Ulix_vs_DMSO_allDEGs.csv")