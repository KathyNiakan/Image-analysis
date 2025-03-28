### Claude

# Load required libraries
library(fgsea)
library(dplyr)
library(ggplot2)
library(tidyr)

# Set work directory
setwd("~/Documents/Collaborations/Brickman/121224")

# # Biololigal process
# pathways <- gmtPathways("~/Documents/Bioinformatics/Databases/m5.go.bp.v2024.1.Mm.symbols.gmt")
# 
# 
# # MF
# pathwaysMF <- gmtPathways("~/Documents/Bioinformatics/Databases/m5.go.mf.v2024.1.Mm.symbols.gmt")

# Hallmark
pathways <- gmtPathways("~/Documents/Bioinformatics/Databases/mh.all.v2024.1.Mm.symbols.gmt")



source("functions/run_enrich.R")



### EPI

# Read DESeq2 results and run GSEA as before
deseq_results <- read.csv("EPI_Ulix_vs_DMSO_allDEGs.csv")
ranks <- deseq_results$stat
names(ranks) <- deseq_results$gene_symbol
ranks <- sort(ranks, decreasing = TRUE)


# Run GSEA
fgsea_results <- fgsea(
  pathways = pathways,
  stats = ranks,
  minSize = 15,
  maxSize = 500,
)

# Filter significant pathways
significant_pathways <- fgsea_results %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>%
  arrange(pval)

# Apply redundancy filtering
filtered_pathways <- filter_redundant_pathways(
  significant_pathways, 
  pathways, 
  similarity_threshold = 0.5  # Adjust this threshold as needed
)

# Format results
results_table <- filtered_pathways %>%
  select(pathway, pval, padj, ES, NES, size) %>%
  mutate(
    pval = format(pval, scientific = TRUE, digits = 2),
    padj = format(padj, scientific = TRUE, digits = 2),
    leadingEdge = sapply(filtered_pathways$leadingEdge, paste, collapse=", ")
  )


# Print summary
cat("Original number of significant pathways:", nrow(significant_pathways), "\n")
cat("Number of non-redundant pathways:", nrow(filtered_pathways), "\n")

# Create visualization of non-redundant pathways
top_pathways <- filtered_pathways %>%
  arrange(desc(abs(NES))) %>%
  head(30)

# Plot
g <- ggplot(top_pathways, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = NES > 0)) +
  coord_flip() +
  labs(x = "Pathway",
       y = "Normalized Enrichment Score",
       title = "Top 30 Non-redundant Enriched Pathways") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(values = c("blue", "red"),
                    name = "Direction",
                    labels = c("Down", "Up"))
print(g)
ggsave("gsea/epi_H.png", g)

# Save filtered results
write.csv(results_table, "gsea/epi_gsea_filtered_results_H.csv", row.names = FALSE)


library(viridis)

g <- ggplot(filtered_pathways, aes(x = NES, y = reorder(pathway, NES), color = padj)) +
  geom_point(size = 4) + scale_color_viridis(option = "magma") +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway",
    title = "Hallmark Pathway Enrichment Analysis"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10)) # Adjust y-axis text size

# Print the plot
print(g)
ggsave("epi_NES_plot.png", width = 8, height = 6, dpi = 300)  # Saves the plot as an image



# Remove "HALLMARK_" prefix from pathway names
filtered_pathways$pathway <- gsub("^HALLMARK_", "", filtered_pathways$pathway)

# Create the bar chart
g <- ggplot(filtered_pathways, aes(x = NES, y = reorder(pathway, NES), fill = padj)) +
  geom_col() +  # Use geom_col() instead of geom_point() for a bar chart
  scale_fill_viridis(option = "magma", name = "Padj", direction = -1) +  # Use fill instead of color
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway",
    title = "Hallmark Pathway Enrichment Analysis"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))  # Adjust y-axis text size

# Print the plot
print(g)
ggsave("epi_NES_plot.png", width = 8, height = 6, dpi = 300)  # Saves the plot as an image



# Create visualization of non-redundant pathways
top_pathways <- filtered_pathways %>%
  arrange(desc(abs(NES))) %>%
  head(5)

# Create the bar chart
g <- ggplot(top_pathways, aes(x = NES, y = reorder(pathway, NES))) +
  geom_col(fill = "purple") +  
  labs(
    x = "NES",
    y = "Pathway"
  ) + theme_classic() + 
  theme(text = element_text(size = 8))  # Adjust y-axis text size

# Print the plot
print(g)
ggsave("gsea/epi_H_top5.png", g, width = 800, height = 300, dpi = 300, units = "px", scale = 1)



### ICM

# Read DESeq2 results and run GSEA as before
deseq_results <- read.csv("ICM_Ulix_vs_DMSO_allDEGs.csv")
ranks <- deseq_results$stat
names(ranks) <- deseq_results$gene_symbol
ranks <- sort(ranks, decreasing = TRUE)


# Run GSEA
fgsea_results <- fgsea(
  pathways = pathways,
  stats = ranks,
  minSize = 5,
  maxSize = 1500,
)

# Filter significant pathways
significant_pathways <- fgsea_results %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>%
  arrange(pval)

# Apply redundancy filtering
filtered_pathways <- filter_redundant_pathways(
  significant_pathways, 
  pathways, 
  similarity_threshold = 0.5  # Adjust this threshold as needed
)

# Format results
results_table <- filtered_pathways %>%
  select(pathway, pval, padj, ES, NES, size) %>%
  mutate(
    pval = format(pval, scientific = TRUE, digits = 2),
    padj = format(padj, scientific = TRUE, digits = 2),
    leadingEdge = sapply(filtered_pathways$leadingEdge, paste, collapse=", ")
  )


# Print summary
cat("Original number of significant pathways:", nrow(significant_pathways), "\n")
cat("Number of non-redundant pathways:", nrow(filtered_pathways), "\n")

# Create visualization of non-redundant pathways
top_pathways <- filtered_pathways %>%
  arrange(desc(abs(NES))) %>%
  head(30)

# Plot
g <- ggplot(top_pathways, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = NES > 0)) +
  coord_flip() +
  labs(x = "Pathway",
       y = "Normalized Enrichment Score",
       title = "Top 30 Non-redundant Enriched Pathways") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(values = c("blue", "red"),
                    name = "Direction",
                    labels = c("Down", "Up"))
print(g)
ggsave("gsea/ICM_H.png", g)

# Save filtered results
write.csv(results_table, "gsea/ICM_gsea_filtered_results_H.csv", row.names = FALSE)
