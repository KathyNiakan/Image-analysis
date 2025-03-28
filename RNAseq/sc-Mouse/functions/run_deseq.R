# To be used to run DESeq2 on cluster by cluster basis (for reannoated clusters)
# Comparison by Treatment: Ulixertinib vs DMSO 
# Run in 02a DEG script

# Function to create and run DESeq2 for a specific cluster

run_deseq_cluster <- function(counts, metadata, cluster_name) {
  # Subset data for specific cluster
  if(cluster_name == "ICM") {
    cluster_cells <- metadata$combined_cluster == cluster_name
  } else {
    cluster_cells <- metadata$cluster_label == cluster_name
  }
  cluster_counts <- counts[, cluster_cells]
  cluster_meta <- metadata[cluster_cells, ]
  
  # Create DESeq dataset
  dds <- DESeqDataSetFromMatrix(
    countData = cluster_counts,
    colData = cluster_meta,
    design = ~Treatment
  )
  
  # Here it would be necessary to filter out no-show genes
  # Already done by Brickman group 
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # VST Transformation
  vsd <- vst(dds, blind = TRUE)
  
  # Extract VST values (normalized counts)
  vsd_mat <- assay(vsd)
  
  # Convert to data frame
  vsd_df <- as.data.frame(vsd_mat)
  vsd_df$gene_id <- rownames(vsd_df)
  
  # Add gene symbols
  vsd_df$gene_symbol <- gene_annotations$gene_symbol[match(vsd_df$gene_id, rownames(gene_annotations))]
  
  # Save VST normalized counts to a CSV file
  vst_file <- paste0("VST_normalized_", cluster_name, ".csv")
  write.csv(vsd_df, vst_file, row.names = FALSE)
  
  # # Save VST normalized counts to a CSV file
  # vst_file <- paste0("VST_normalized_", cluster_name, ".csv")
  # write.csv(vsd_mat, vst_file)
  
  # Get results (Ulixertinib vs DMSO)
  res <- results(dds, contrast = c("Treatment", "Ulixertinib", "DMSO"))
  
  # Convert to data frame and add gene info
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df$cluster <- cluster_name
  
  return(res_df)
}


# Get significant DEGs for each cluster
get_sig_degs <- function(deg_df, padj_threshold = 0.05, lfc_threshold = 1.5) {
  subset(deg_df, padj < padj_threshold & abs(log2FoldChange) > lfc_threshold)
}

# Get all DEGs for each cluster for volcano
get_all_degs <- function(deg_df) {
  subset(deg_df, padj > 0 & abs(log2FoldChange) > 0)
}

