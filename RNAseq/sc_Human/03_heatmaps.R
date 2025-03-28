# Heatmaps ####

#===============================================================================
# Load packages, set wd
#===============================================================================


# Required packages
library(pheatmap)
library(dplyr)
library(tidyr)
library(readr)
library(viridis)

# Set work directory
setwd("~/Documents/Collaborations/Brickman/61224/11_niakan_adata_combined")


#===============================================================================
# EPI Cells load data and wrangling
#===============================================================================

# Read in data
EPI_counts <- read.csv("VST_normalized_EPI.csv")

info <- c("gene_id","gene_symbol")
gene_info <- EPI_counts[, info]

# Get gene annotations and ensure symbols are character, not factor
symbols <- as.character(gene_info$gene_symbol)

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
rownames(EPI_counts) <- unique_symbols

# Drop extra columns
EPI_counts <- EPI_counts[, !(colnames(EPI_counts) %in% info)]

# Convert counts to a numeric matrix
expression_matrix <- as.matrix(EPI_counts)


metadata <- read.csv("cell_metadata_cs.csv") 
rownames(metadata) <- metadata[, 1]
metadata <- metadata[, -1]


# Check alignment between metadata and counts
if (!all(colnames(expression_matrix) %in% rownames(metadata))) {
  stop("Mismatch between cells in counts and metadata.")
}

# Subset to matching cells
expression_matrix <- expression_matrix[, colnames(expression_matrix) %in% rownames(metadata)]
metadata <- metadata[rownames(metadata) %in% colnames(expression_matrix), ]



#===============================================================================
# EPI Cells Heatmap
#===============================================================================

# Function
create_expression_heatmap <- function(expression_matrix, metadata, genes_of_interest) {
  # Subset expression matrix to genes of interest
  subset_matrix <- expression_matrix[genes_of_interest, ]
  
  # Remove genes with zero variance
  row_vars <- apply(subset_matrix, 1, var, na.rm=TRUE)
  subset_matrix <- subset_matrix[row_vars > 0, ]
  
  # If any genes were removed, print a warning
  removed_genes <- genes_of_interest[row_vars == 0]
  if(length(removed_genes) > 0) {
    warning("Removed genes with zero variance: ", paste(removed_genes, collapse=", "))
  }
  
  # Create annotation dataframe for treatments
  annotation_df <- data.frame(
    Treatment = metadata$treatment,
    Cluster = metadata$cluster_label,
    row.names = rownames(metadata)
  )
  
  annotation_colors <- list(
    Treatment = c(
      "Ulix" = "orange",
      "DMSO" = "black"),
    Cluster = c("EPI" = "#00FF00")
  )
  
  breaks <- seq(-2, 2, length.out = 10)
  
  pheatmap(
    subset_matrix,
    annotation_col = annotation_df,
    annotation_colors = annotation_colors,
    show_rownames = TRUE,
    show_colnames = FALSE,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    treeheight_row = 0,      
    treeheight_col = 5,  
    scale = "row",
    height = 3,
    width = 3.5,
    fontsize = 8,
    color = inferno(10),
    breaks = breaks,
    filename = "EPI_heatmap.png" # also .pdf for changing to italics - can't do this in Pheatmap
  )
}

epi_deg_OI <- (read.csv("EPI_genelist.csv", header = FALSE))
gene_list <- epi_deg_OI[,1]

# Usage example:
create_expression_heatmap(
  expression_matrix = expression_matrix,
  metadata = metadata,
  genes_of_interest = gene_list
)



#===============================================================================
# PrE Cells load data and wrangling
#===============================================================================


# Read in data
PrE_counts <- read.csv("VST_normalized_PrE.csv")

info <- c("gene_id","gene_symbol")
gene_info <- PrE_counts[, info]

# Get gene annotations and ensure symbols are character, not factor
symbols <- as.character(gene_info$gene_symbol)

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
rownames(PrE_counts) <- unique_symbols

# Drop extra columns
PrE_counts <- PrE_counts[, !(colnames(PrE_counts) %in% info)]

# Convert counts to a numeric matrix
expression_matrix <- as.matrix(PrE_counts)


metadata <- read.csv("cell_metadata_cs.csv") 
rownames(metadata) <- metadata[, 1]
metadata <- metadata[, -1]


# Check alignment between metadata and counts
if (!all(colnames(expression_matrix) %in% rownames(metadata))) {
  stop("Mismatch between cells in counts and metadata.")
}

# Subset to matching cells
expression_matrix <- expression_matrix[, colnames(expression_matrix) %in% rownames(metadata)]
metadata <- metadata[rownames(metadata) %in% colnames(expression_matrix), ]



#===============================================================================
# PrE Cells Heatmap
#===============================================================================

# Function
create_expression_heatmap <- function(expression_matrix, metadata, genes_of_interest) {
  # Subset expression matrix to genes of interest
  subset_matrix <- expression_matrix[genes_of_interest, ]
  
  # Remove genes with zero variance
  row_vars <- apply(subset_matrix, 1, var, na.rm=TRUE)
  subset_matrix <- subset_matrix[row_vars > 0, ]
  
  # If any genes were removed, print a warning
  removed_genes <- genes_of_interest[row_vars == 0]
  if(length(removed_genes) > 0) {
    warning("Removed genes with zero variance: ", paste(removed_genes, collapse=", "))
  }
  
  # Create annotation dataframe for treatments
  annotation_df <- data.frame(
    Treatment = metadata$treatment,
    Cluster = metadata$cluster_label,
    row.names = rownames(metadata)
  )
  
  annotation_colors <- list(
    Treatment = c(
      "Ulix" = "orange",
      "DMSO" = "black"),
    Cluster = c("PrE" = "magenta")
  )
  
  breaks <- seq(-2, 2, length.out = 10)
  
  pheatmap(
    subset_matrix,
    annotation_col = annotation_df,
    annotation_colors = annotation_colors,
    show_rownames = TRUE,
    show_colnames = FALSE,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    treeheight_row = 0,      
    treeheight_col = 5,  
    scale = "row",
    height = 3,
    width = 3.5,
    fontsize = 8,
    color = inferno(10),
    breaks = breaks,
    filename = "PrE_heatmap.pdf"
  )
}

PrE_deg_OI <- (read.csv("PrE_genelist.csv", header = FALSE))
gene_list <- PrE_deg_OI[,1]

# Usage example:
create_expression_heatmap(
  expression_matrix = expression_matrix,
  metadata = metadata,
  genes_of_interest = gene_list
)







#===============================================================================
# TE Cells load data and wrangling
#===============================================================================

# Read in data
TE_counts <- read.csv("VST_normalized_TE.csv")

metadata <- read.csv("cell_metadata_cs.csv") 

TE_counts <- TE_counts[!duplicated(TE_counts[, 1]), ]
rownames(TE_counts) <- TE_counts[, 1]
TE_counts <- TE_counts[, -1]

rownames(metadata) <- metadata[, 1]
metadata <- metadata[, -1]

# Convert counts to a numeric matrix
expression_matrix <- as.matrix(TE_counts)

# Check alignment between metadata and counts
if (!all(colnames(expression_matrix) %in% rownames(metadata))) {
  stop("Mismatch between cells in counts and metadata.")
}

# Subset to matching cells
expression_matrix <- expression_matrix[, colnames(expression_matrix) %in% rownames(metadata)]
metadata <- metadata[rownames(metadata) %in% colnames(expression_matrix), ]

#===============================================================================
# TE Cells Heatmap
#===============================================================================

# Function
create_expression_heatmap <- function(expression_matrix, metadata, genes_of_interest) {
  # Subset expression matrix to genes of interest
  subset_matrix <- expression_matrix[genes_of_interest, ]
  
  # Remove genes with zero variance
  row_vars <- apply(subset_matrix, 1, var, na.rm=TRUE)
  subset_matrix <- subset_matrix[row_vars > 0, ]
  
  # If any genes were removed, print a warning
  removed_genes <- genes_of_interest[row_vars == 0]
  if(length(removed_genes) > 0) {
    warning("Removed genes with zero variance: ", paste(removed_genes, collapse=", "))
  }
  
  # Create annotation dataframe for treatments
  annotation_df <- data.frame(
    Treatment = metadata$treatment,
    Cluster = metadata$cluster_label,
    row.names = rownames(metadata)
  )
  
  annotation_colors <- list(
    Treatment = c(
      "Ulix" = "orange",
      "DMSO" = "black"),
    Cluster = c("TE" = "#0000FF")
  )
  
  breaks <- seq(-2, 2, length.out = 10)
  
  pheatmap(
    subset_matrix,
    annotation_col = annotation_df,
    annotation_colors = annotation_colors,
    show_rownames = TRUE,
    show_colnames = FALSE,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    treeheight_row = 0,      
    treeheight_col = 5,  
    scale = "row",
    height = 2,
    width = 4,
    fontsize = 5,
    color = inferno(10),
    breaks = breaks,
    filename = "TE_heatmap.png"
  )
}

TE_deg_OI <- (read.csv("TE_genelist.csv", header = FALSE))
gene_list <- TE_deg_OI[,1]

# Usage example:
create_expression_heatmap(
  expression_matrix = expression_matrix,
  metadata = metadata,
  genes_of_interest = gene_list
)









#===============================================================================
# ICM Cells load data and wrangling
#===============================================================================

# Read in data
ICM_counts <- read.csv("VST_normalized_ICM.csv")

metadata <- read.csv("cell_metadata_cs.csv") 

ICM_counts <- ICM_counts[!duplicated(ICM_counts[, 1]), ]
rownames(ICM_counts) <- ICM_counts[, 1]
ICM_counts <- ICM_counts[, -1]

rownames(metadata) <- metadata[, 1]
metadata <- metadata[, -1]

# Convert counts to a numeric matrix
expression_matrix <- as.matrix(ICM_counts)

# Check alignment between metadata and counts
if (!all(colnames(expression_matrix) %in% rownames(metadata))) {
  stop("Mismatch between cells in counts and metadata.")
}

# Subset to matching cells
expression_matrix <- expression_matrix[, colnames(expression_matrix) %in% rownames(metadata)]
metadata <- metadata[rownames(metadata) %in% colnames(expression_matrix), ]

#===============================================================================
# ICM Cells Heatmap
#===============================================================================

# Function
create_expression_heatmap <- function(expression_matrix, metadata, genes_of_interest) {
  # Subset expression matrix to genes of interest
  subset_matrix <- expression_matrix[genes_of_interest, ]
  
  # Remove genes with zero variance
  row_vars <- apply(subset_matrix, 1, var, na.rm=TRUE)
  subset_matrix <- subset_matrix[row_vars > 0, ]
  
  # If any genes were removed, print a warning
  removed_genes <- genes_of_interest[row_vars == 0]
  if(length(removed_genes) > 0) {
    warning("Removed genes with zero variance: ", paste(removed_genes, collapse=", "))
  }
  
  # Create annotation dataframe for treatments
  annotation_df <- data.frame(
    Treatment = metadata$treatment,
    Cluster = metadata$cluster_label,
    row.names = rownames(metadata)
  )
  
  annotation_colors <- list(
    Treatment = c(
      "Ulix" = "orange",
      "DMSO" = "black"), 
    Cluster = c("EPI" = "#00FF00",
                                     "PrE" = "magenta")
  )
  
  breaks <- seq(-2, 2, length.out = 10)
  
  pheatmap(
    subset_matrix,
    annotation_col = annotation_df,
    annotation_colors = annotation_colors,
    show_rownames = TRUE,
    show_colnames = FALSE,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    treeheight_row = 0,      
    treeheight_col = 5,  
    scale = "row",
    height = 3,
    width = 4,
    fontsize = 5,
    color = inferno(10),
    breaks = breaks,
    filename = "ICM_heatmap.png"
  )
}

ICM_deg_OI <- (read.csv("ICM_genelist.csv", header = FALSE))
gene_list <- ICM_deg_OI[,1]

# Usage example:
create_expression_heatmap(
  expression_matrix = expression_matrix,
  metadata = metadata,
  genes_of_interest = gene_list
)
