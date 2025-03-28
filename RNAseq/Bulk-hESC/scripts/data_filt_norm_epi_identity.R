# Data filtering and normalisation
# Deconvolution of lineage identity 
# Plot results Fig 4H

setwd("/Users/simonc/Documents/Cell_work/Human/RNAseq/scripts")

library(RColorBrewer)
library(Rtsne)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(readr)
library(DrImpute)
library(purrr)
library(tidyr)
library(DeconRNASeq)
library(scales)
library(stringr)
library(tidyverse)
library(tximport)

source("utility_functions.R")

## Skip to line 257 to read in processed results and re-plot

# Load sample metadata
stype <- read_csv("../metadata/more_samples_features.csv")
stype <- stype[!duplicated(stype$sample_name),]
# Fix weird characters in cell_line
stype <- stype %>% mutate(cell_line = if_else(cell_line == "na\x95ve_hESC", "naive_hESC", cell_line))
#stype <- stype %>% mutate(sample_name = str_replace(sample_name, "-", "."))

# Load the hESC and additional study count data
load("../results/hesc.RData")
names(hesc)
load("../results/hesc_new.RData")
names(hesc_new)

# Function to combine technical replicates by summing (for counts)
combine_technical_replicates_sum <- function(count_matrix) {
  summed_matrix <- count_matrix %>%
    as.data.frame() %>%
    rownames_to_column("Gene") %>%
    pivot_longer(-Gene, names_to = "Sample", values_to = "Count") %>%
    mutate(Condition = sub("-rep[0-9]+$", "", Sample)) %>%
    group_by(Gene, Condition) %>%
    summarise(Summed_Count = sum(Count), .groups = "drop") %>%
    pivot_wider(names_from = Condition, values_from = Summed_Count)
  
  # Convert back to a matrix with row names
  summed_matrix <- summed_matrix %>%
    column_to_rownames("Gene")
  
  return(as.matrix(summed_matrix))
}

# Function to combine technical replicates by averaging (for TPM and lengths)
combine_technical_replicates_mean <- function(matrix_data) {
  averaged_matrix <- matrix_data %>%
    as.data.frame() %>%
    rownames_to_column("Gene") %>%
    pivot_longer(-Gene, names_to = "Sample", values_to = "Value") %>%
    mutate(Condition = sub("-rep[0-9]+$", "", Sample)) %>%
    group_by(Gene, Condition) %>%
    summarise(Avg_Value = mean(Value), .groups = "drop") %>%
    pivot_wider(names_from = Condition, values_from = Avg_Value)
  
  # Convert back to a matrix with row names
  averaged_matrix <- averaged_matrix %>%
    column_to_rownames("Gene")
  
  return(as.matrix(averaged_matrix))
}

# Process each component of the hesc_new object
hesc_new_processed <- list()

# Process counts (sum)
hesc_new_processed$counts <- combine_technical_replicates_sum(hesc_new$counts)

# Process abundance (TPM values) - using mean instead of sum
hesc_new_processed$abundance <- combine_technical_replicates_mean(hesc_new$abundance)

# Process length (using mean)
hesc_new_processed$length <- combine_technical_replicates_mean(hesc_new$length)

# Set countsFromAbundance to match original setting
hesc_new_processed$countsFromAbundance <- hesc_new$countsFromAbundance

# Save the processed data
save(hesc_new_processed, file = "../results/hesc_new_processed.RData")

hesc_new <- hesc_new_processed

# Subset and save claire counts and tpm for geo submission
geo_counts <- hesc$counts
geo_counts <- geo_counts[,grep("SIM", colnames(geo_counts))]

hesc_new_counts <- hesc_new$counts 

# 
# # Fix metadata to reflect summed
# stype_hesc_new <- stype_hesc_new %>%
#   mutate(sample_name = sub("-rep[0-9]+$", "", sample_name)) %>%
#   distinct(sample_name, .keep_all = TRUE)

geo_counts <- merge(geo_counts, hesc_new_counts, by = "row.names")
rownames(geo_counts) <- geo_counts$Row.names
geo_counts$Row.names <- NULL

geo_counts <- data.frame(geo_counts)
geo_counts$gene <- rownames(geo_counts)
geo_counts <- geo_counts[,c(11,1:10)]


write_tsv(geo_counts, "../results/geo_counts.tsv")

geo_tpm <- hesc$abundance
geo_tpm <- geo_tpm[,grep("SIM", colnames(geo_tpm))]
hesc_new_tpm <- hesc_new$abundance

geo_tpm <- merge(geo_tpm, hesc_new_tpm, by = "row.names")
rownames(geo_tpm) <- geo_tpm$Row.names
geo_tpm$Row.names <- NULL

geo_tpm <- data.frame(geo_tpm)
geo_tpm$gene <- rownames(geo_tpm)
geo_tpm <- geo_tpm[,c(11,1:10)]
write_tsv(geo_tpm, "../results/geo_tpm.tsv")

# Load Wamaitha datasets
wamaitha <- readRDS("../results/wamaitha.RDS")

# Extract count/length data
wam_counts <- assays(wamaitha)$counts
wam_len <- assays(wamaitha)$avgTxLength

# Calculate TPM from counts data
tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}
hesc_tpm <- tpm3(hesc$counts, hesc$length)
hesc_new_tpm <- tpm3(hesc_new$counts, hesc_new$length)
wam_tpm <- tpm3(wam_counts, wam_len)

# Merge together tpm tables
tpm_allhESC <- merge(hesc_tpm, hesc_new_tpm, by = "row.names")
rownames(tpm_allhESC) <- tpm_allhESC$Row.names
tpm_allhESC <- tpm_allhESC[,-1]
tpm <- merge(tpm_allhESC, wam_tpm, by = "row.names")

rownames(tpm) <- tpm$Row.names
tpm$Row.names <- NULL

# Limit metadata to samples present in tpm data ### here is where the samples in stype go missing
stype <- stype[stype$sample_name %in% colnames(tpm),]


tpm <- tpm[,stype$sample_name]
all(colnames(tpm) == stype$sample_name)


# Now get single-cell epi, te, pe data
# Load sample metadata
meta <- read_tsv("../metadata/sample_features.tsv")
meta_new <- read_csv("../metadata/new_sample_features.csv")

meta <- rbind(meta, meta_new)

meta <- meta[!duplicated(meta$sample_name),]
# Fix weird characters in cell_line
meta <- meta %>% mutate(cell_line = if_else(cell_line == "na\x95ve_hESC", "naive_hESC", cell_line))



# Limit metadata to samples present in hesc counts
meta_ref <- meta[meta$cell_line %in% c("Epi", "TE", "PE"),]
# 
# 
# wamaitha <- readRDS("results/wamaitha.RDS") # Already in summarized experiment format
# 
# counts_wam <- assays(wamaitha)$counts
# len_wam <- assays(wamaitha)$avgTxLength
# tpm_wam <- tpm3(counts_wam, len_wam)
# stype_wam <- meta_ref[meta_ref$sample_name %in% colnames(tpm_wam),]
# 
# table(stype_wam$batch)

ref_tpm <- wam_tpm[,meta_ref$sample_name]
all(colnames(ref_tpm) == meta_ref$sample_name)

# tpm <- merge(tpm_wam, hesc_tpm, by = "row.names")
# rownames(tpm) <- tpm$Row.names
# tpm$Row.names <- NULL
# 
# ref_tpm <- ref_tpm[ref_tpm$sample_name %in% colnames(tpm_wam),]


# Split out epi, te, pre data to generate reference expression set 
meta_ref_epi <- meta_ref[meta_ref$cell_line == "Epi",]
meta_ref_te <- meta_ref[meta_ref$cell_line == "TE",]
meta_ref_pe <- meta_ref[meta_ref$cell_line == "PE",]

tpm_epi <- ref_tpm[,meta_ref_epi$sample_name]
tpm_te <- ref_tpm[,meta_ref_te$sample_name]
tpm_pe <- ref_tpm[,meta_ref_pe$sample_name]

epi_mean <- data.frame(gene = rownames(tpm_epi), epi = rowMeans(tpm_epi))
te_mean <- data.frame(gene = rownames(tpm_te), te = rowMeans(tpm_te))
pe_mean <- data.frame(gene = rownames(tpm_pe), pe = rowMeans(tpm_pe))

tpm_list <- list(epi_mean, te_mean, pe_mean)

tpm_blast <- tpm_list %>% reduce(merge, by='gene')
rownames(tpm_blast) <- tpm_blast$gene
tpm_blast$gene <- NULL

# Remove genes with tpm < 10 in any of the reference tissues
# Results in 5494 genes
tpm_blast_filt <- tpm_blast[!apply(tpm_blast<10,1,any),]


# Get tpm for test samples - Claire hESC
# claire_stype <- stype[stype$batch == "Simon",]
# claire_tpm <- tpm[,claire_stype$sample_name]

# Run deconvolution
signatures <- tpm_blast_filt
datasets <- tpm

DeconRNASeq(datasets, signatures, checksig=FALSE,
            known.prop = F, use.scale = TRUE, fig = TRUE)
res <- DeconRNASeq(datasets, signatures, checksig=TRUE, 
                   known.prop = F, use.scale = TRUE, fig = TRUE)


# Format deconvolution results
res_prop <- res$out.all
res_prop <- data.frame(res_prop)
res_prop$sample <- stype$medium
res_prop$state <- stype$state

res_prop_long <- res_prop %>% pivot_longer(c(epi, te, pe), names_to = "cell_type", values_to = "proportion")

res_prop_long <- res_prop_long %>% 
  mutate(treatment = case_when(res_prop_long$sample == "DMSO" ~ "Control in PXGL",
                               res_prop_long$sample == "ULIX" ~ "ERKi in PXGL",
                               res_prop_long$sample == "UXGL_DMSO" ~ "Control in UXGL",
                               res_prop_long$sample == "UXGL_ULIX" ~ "ERKi in UXGL",
                               TRUE ~ state))


write.csv(res_prop_long, "../results/percentage_epi_data.csv", row.names = F)

# Read in res_prop_long data instead of repeating upstream
res_prop_long <- read.csv("../results/percentage_epi_data.csv")

res_prop_long$treatment <- factor(res_prop_long$treatment, levels = c("primed", "naive", "Control in PXGL", "ERKi in PXGL", "Control in UXGL", "ERKi in UXGL"),
                                  labels = c("Primed hESC", "Naive hESC", "Control in PXGL", "ERKi in PXGL", "Control in UXGL", "ERKi in UXGL"))


# Plot results
p <- ggplot(res_prop_long[res_prop_long$cell_type == "epi",], aes(x = treatment, y = proportion, fill = treatment)) + 
  theme_classic(base_size = 8)
p <- p + stat_summary (aes (color = treatment), position = position_dodge(width = 0.7), fun.y = median, fun.ymin = median, 
                         fun.ymax = median, geom = "crossbar", width = 0.7, size = 0.7)
p <- p + geom_jitter (aes (color = treatment), 
                     position = position_dodge(width = 0.7), alpha = 0.3,  size = 1, 
                     stroke = 0.5) 
p <- p +   scale_color_manual(values = c("#00B865", "#70A7FF", "#9C9996", "#EF850B", "#000000", "#FB4570")) + scale_fill_manual(values = c("#00B865", "#70A7FF", "#9C9996", "#EF850B",   "#000000", "#FB4570"))  + theme_classic() +
  scale_fill_manual(values = c("#00B865", "#70A7FF", "#9C9996", "#EF850B","#000000", "#FB4570")) +
  labs(y = "% Epi-like \n identity", x = "") + 
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(guide = guide_axis(angle = 45))+theme(aspect.ratio = 0.6, legend.position="none")
print(p)
ggsave("../results/percentage_epi_boxplot.png", p, height = 2.2, width = 3)
ggsave("../results/percentage_epi_boxplot.pdf", p, height = 3)

