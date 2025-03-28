# Heatmaps ####
# ==============================================================================
# Load packages, set wd
# ==============================================================================

# Required packages
library(pheatmap)
library(dplyr)
library(tidyr)
library(readr)
library(viridis)


# Set work directory
setwd("~/Documents/Collaborations/Brickman/121224")

# Load data
load("sce_filtered_cscluster.RData")

# Get gene annotations and ensure symbols are character, not factor
gene_annotations <- as.data.frame(rowData(sce))
symbols <- as.character(gene_annotations$gene_symbol)

# Extract the counts matrix
expr <- assay(sce, "logcounts")
expr_df <- as.data.frame(as.matrix(expr))

# Find duplicated symbols
dups <- duplicated(symbols) | duplicated(symbols, fromLast = TRUE)
dup_symbols <- unique(symbols[dups])

# Create unique names
unique_symbols <- symbols
for(sym in dup_symbols) {
  # Find all instances of this symbol
  idx <- which(symbols == sym)
  # Replace with numbered versions
  unique_symbols[idx] <- paste0(sym, "-", seq_along(idx))
}

# Assign as rownames
rownames(expr_df) <- unique_symbols

# Convert back to matrix
expr <- as.matrix(expr_df)


metadata <- as.data.frame(colData(sce))
metadata$sample_name <- rownames(metadata)


# Check alignment between metadata and counts
if (!all(colnames(expr) %in% rownames(metadata))) {
  stop("Mismatch between cells in RPM counts and metadata.")
}

# Subset to matching cells
expr <- expr[, colnames(expr) %in% rownames(metadata)]
metadata <- metadata[rownames(metadata) %in% colnames(expr), ]

# ==============================================================================
# Heatmap for all cells
# ==============================================================================

# 1. Basic heatmap equivalent to the clustermap
create_expression_heatmap <- function(expression_matrix, treatment_metadata, genes_of_interest) {
  # Subset expression matrix to genes of interest
  subset_matrix <- expression_matrix[genes_of_interest, ]
  
  # Create annotation dataframe for treatments
  annotation_df <- data.frame(
    Treatment = treatment_metadata$Treatment,
    Brickman = treatment_metadata$predictions,
    Simon = treatment_metadata$cluster_label,
    row.names = rownames(treatment_metadata)
  )
  
  # Create color palette for treatments
  annotation_colors <- list(
    Treatment = c(
      "Ulixertinib" = "orange",
       "DMSO" = "black"),
    Brickman = c("E3.5-ICM" = "#7F7F7F",
                 "E4.5-EPI" = "#00FF00", 
                 "E4.5-PrE" = "#FF00FF"),
    Simon = c("EPI" = "#00FF00",
               "PrE" = "#FF00FF")
  )
#  breaks <- seq(-2, 2, length.out = 10)  # Adjust breaks for scaling, modify the range (-2 to 2) as needed
  
  # Create heatmap
  pheatmap(
    subset_matrix,
    annotation_col = annotation_df,
    annotation_colors = annotation_colors,
    show_rownames = TRUE,
    show_colnames = FALSE,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
   # scale = "row",  # Scale by row (gene)
    height = 8,
    width = 5,
    color = inferno(10),
  # breaks = breaks
    
  )
}

gene_list <- c("Nanog", "Sox2", "Klf4", "Fgf4", "Etv5",
              "Pou5f1", "Dab2",
               "Pdgfra", "Gata6", "Fgfr2",
               "Gata4",  "Sox17", "Sox7")

# EPI_goi <- read.csv("EPI_goi.csv", header = FALSE)
# gene_list <- EPI_goi$V1


# Usage example:
create_expression_heatmap(
  expression_matrix = expr,
  treatment_metadata = metadata,
  genes_of_interest = gene_list
)
