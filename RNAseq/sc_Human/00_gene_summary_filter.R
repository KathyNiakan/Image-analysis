# This is code to use with raw data counts of human single-cell SMART-seq
# Code is based on https://github.com/brickmanlab/proks-salehin-et-al and 
# https://github.com/galanisl/AI_hESCs and from Laura Woods


# ==============================================================================
# Load packages, set wd, load data
# ==============================================================================
library(biomaRt)
library(SingleCellExperiment)
library(scran)
library(dplyr)
library(stringr)
library(zellkonverter)
library(readr)


# Set work directory
setwd("~/Documents/Collaborations/Brickman/61224/11_niakan_adata_combined")

# Raw counts, not length normalised
sce <- readH5AD("11_niakan_adata_combined_notLengthNormalised_withPredictions.h5ad")

# Rename X (Python) to counts (R)
assayNames(sce) <- c("counts")

# ==============================================================================
# Calculate TPM and FPKM
# ==============================================================================

# Get gene lengths for Ensembl gene IDs via BioMart
# Takes a couple of minutes
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_ids <- rownames(sce)

gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "transcript_length", "chromosome_name"
                                  # , 
                                  # "gene_biotype" # Include biotype for filtering pseudo genes
                                  ),
                   filters = "ensembl_gene_id",
                   values = ensembl_ids,
                   mart = mart) %>%
  group_by(ensembl_gene_id) %>%
  summarize(
    gene_name = first(external_gene_name),  
    length = median(transcript_length),      
    chr = first(chromosome_name),
#    gene_biotype = first(gene_biotype)
  ) %>%
  as.data.frame()

gene_info$gene_name[gene_info$gene_name == ""] <- NA
gene_info$gene_name[is.na(gene_info$gene_name)] <- gene_info$ensembl_gene_id[is.na(gene_info$gene_name)]

gene_info %>% filter(length == "")


# TPM-like

counts_matrix <- counts(sce)
libSize <- colSums(counts_matrix)

# Function to perform gene length normalization
gene_lengths <- gene_info$length[match(rownames(counts_matrix), gene_info$ensembl_gene_id)]
rpk <- counts_matrix / (gene_lengths/1000)

# Main normalization workflow for SCE object
scaling_factors <- libSize / 1e6
tpm_like <- t(t(rpk) / scaling_factors)
log_normalized <- log1p(tpm_like)
  
# Store all normalizations as assays
assays(sce)$length_normalized <- rpk
assays(sce)$tpm <- tpm_like
assays(sce)$logtpm <- log_normalized


# Calculate FPKM
fpkm <- (counts_matrix / (gene_lengths * matrix(libSize, 
                                                   byrow = TRUE, 
                                                   nrow(counts_matrix), 
                                                   ncol(counts_matrix)))) * 1e9
# Combine to SCE
assays(sce)$fpkm <- fpkm

# ==============================================================================
# Gene-level filtering
# ==============================================================================

# Note that sample-level filtering of lowly expressing samples and 
# high Mt content already done by Brickman lab

# Get ensembl ids and gene names
genes_sce <- as.data.frame(rowData(sce))

# filter ribosomal genes
ribo <- read_tsv("ribosomal_human.txt")
ribo_sce <- genes_sce %>% filter(gene_symbol %in% ribo$`Approved Symbol`)
ribo_drop <- rownames(ribo_sce)
sce <- sce[!(rownames(sce) %in% ribo_drop), ]

# remove mitochondrial genes
mito <- gene_info %>% filter(chr == "MT")
sce <- sce[!(rownames(sce) %in% mito$ensembl_gene_id), ]

# Get rid of no-show genes
sce <- sce[rowSums(counts(sce)) != 0, ]

# Only keep genes with avg. expression across cells >= 1
ave.counts <- rowMeans(counts(sce))
keep <- ave.counts >= 1
sum(keep)
sce <- sce[keep, ]


# ==============================================================================
# save file
# ==============================================================================

save(sce, file = "sce_filtered.RData")
