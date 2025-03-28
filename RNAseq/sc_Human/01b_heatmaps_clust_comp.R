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
setwd("~/Documents/Collaborations/Brickman/61224/11_niakan_adata_combined")

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
lanner <- read.table("Lanner/EEB.2024-12-11.anno.tsv", header = TRUE)
colnames(lanner)[colnames(lanner) == "query_cell"] <- "sample_name"
lanner <- lanner[,c(1,2)]

metadata <- left_join(metadata, lanner, by="sample_name")
rownames(metadata) <- metadata$sample_name

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
    Treatment = treatment_metadata$treatment,
    Brickman = treatment_metadata$prediction,
    Lanner = treatment_metadata$pred_EML,
    Simon = treatment_metadata$cluster_label,
    row.names = rownames(treatment_metadata)
  )
  
  # Create color palette for treatments
  annotation_colors <- list(
    Treatment = c(
      "Ulix" = "orange",
       "DMSO" = "black"),
    Brickman = c("Inner Cell Mass" = "#7F7F7F",
                "Epiblast_6.0" = "#00FF00",
                "Epiblast_7.0" = "#00BF00",
                "Late epiblast" = "#007F00",
                "Primitive Endoderm" = "#FF00FF",
                "Trophectoderm_6.0" = "#0000FF",
                "Trophectoderm_7.0" = "#0000BF",
                "Trophectoderm_8.0" = "#00007F",
                "Trophectoderm_9.0" = "#00003F"),
    Lanner = c("ICM" = "#7F7F7F",
                 "Epiblast" = "#00FF00",
                 "PriS" = "#007F00",
                 "Hypoblast" = "#FF00FF",
                 "Ambiguous" = "#BA00BF",
                 "TE" = "#0000FF"),
    Simon = c("EPI" = "#00FF00",
               "PrE" = "#FF00FF",
               "TE" = "#0000FF")
  )
  breaks <- seq(-1, 2.5, length.out = 10)  # Adjust breaks for scaling, modify the range (-2 to 2) as needed
  
  # Create heatmap
  pheatmap(
    subset_matrix,
    annotation_col = annotation_df,
    annotation_colors = annotation_colors,
    show_rownames = TRUE,
    show_colnames = FALSE,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    scale = "row",  # Scale by row (gene)
    height = 8,
    width = 5,
    color = inferno(10),
    breaks = breaks
    
  )
}

gene_list <- c("NANOG", "SOX2", "KLF17", "TDGF1", "PDGFRA", "GATA6", 
               "GATA4", "SOX17", "GATA2", "GATA3", "KRT18", "TEAD3")



# Usage example:
create_expression_heatmap(
  expression_matrix = expr,
  treatment_metadata = metadata,
  genes_of_interest = gene_list
)
