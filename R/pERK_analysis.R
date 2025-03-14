#------------------------------------------------------------------------------#
#  This script accompanies the manuscript                                      #
#  Simon et al., (2025)                                                        #
#  Repository available on:                                                    # 
#  https://ctr.uniofcam.dev/ctr-bioinformatics/niakan-lab/simon-et-al-2025     #
#  Please consult READ_ME for more information                                 #
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#                                pERK: Analysis                                #
#------------------------------------------------------------------------------#

#..............................................................................#
# Input:    Embryos_clust from pERK_data.R                                     #
#                                                                              #
# Outputs:  Fig. 2b                                                            #
#           Fig. 2c                                                            #
#..............................................................................#

# Load library packages, data and set wd ----

library(dplyr)
library(ggplot2)
library(tidyverse)
library(viridis)
library(RColorBrewer)

setwd("/Users/simonc/Documents/Human_Embryos/Simon_et_al_2025/R")

source("pERK_data.R")

# Data wrangling ----

# Add staging information 
Embryos_clust$Stage = ifelse(grepl("220120", Embryos_clust$Embryo_ID), "Day 5",
                             ifelse(grepl("211021", Embryos_clust$Embryo_ID), "Day 6",
                                    "Day 6.5"))

# Count number of cells in each lineage per embryo
Embryos_count <- Embryos_clust %>%
  group_by(Embryo_ID, id.cluster) %>%
  summarize(N = n()) %>%
  ungroup() %>%
  complete(Embryo_ID, id.cluster,
           fill = list(N = 0))

# Total cell counts per embryo
CellNo <- Embryos_clust %>%
  group_by(Embryo_ID) %>%
  summarize(CellNumber = n())

Embryos_count <- Embryos_count %>% left_join(CellNo, by = "Embryo_ID")

# ICM Clustering ----

# Scatter plot of SOX2 and OTX2 nuclear fluorescence intensity
# Check hierarchical clustering assignment to epiblast and hypoblast
g <- Embryos_clust %>% filter(TE_ICM == "ICM")  %>%
  ggplot(aes(x = SOX2_Cor, y = OTX2_Cor))
g <- g + geom_point(aes(fill=id.cluster), 
                    size = 2, shape = 21, colour = "black") 
g <- g +  ylab("log[OTX2]") +  xlab("log[SOX2]") + theme_classic()  + 
  theme(text = element_text(size = 24)) + scale_fill_manual(values = col)
print(g)

# Scatter plot of SOX2 and OTX2 nuclear fluorescence intensity
# Overlay cytoplasmic pERK intensity
g <- Embryos_clust %>% filter(TE_ICM == "ICM")  %>%
  ggplot(aes(x = SOX2_Cor, y = OTX2_Cor))
g <- g + geom_point(aes(fill = Cyto_pERK_Cor), 
                    size = 2, shape = 21, colour = "black") 
g <- g +  ylab("log[OTX2]") +  xlab("log[SOX2]") + theme_classic() + 
  scale_fill_viridis()
print(g)

# Colour palettes ----
col <- c("EPI" = "red", "PrE" = "cyan", "TE" = "grey")

# Fig. 2b: pERK by stage ----

# Violin and boxplot of cytoplasmic pERK intensity, by stage
g <- Embryos_clust %>% 
  ggplot(aes(x = Stage, y = Cyto_pERK_Cor, fill = Stage))
g <- g + geom_violin(trim = FALSE, alpha = 0.7) 
g <- g + geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.1) 
g <- g + theme_classic() +
  ylab("Cytoplasmic log[pERK]") + 
  theme(text = element_text(size = 24),
        axis.text.x = element_text(angle = 90)) + 
  scale_fill_brewer(palette = "BuPu") + scale_color_brewer(palette = "BuPu")
print(g)

ggsave(filename = "figures/Fig2b_pERK_stage.png", g,
       width = 450, height = 450, dpi = 300, units = "px", scale = 4)

# Fig. 2c: pERK by lineage ----

# Violin and boxplot of cytoplasmic pERK intensity, by lineage
g <- Embryos_clust %>% 
  ggplot(aes(x = id.cluster, y = Cyto_pERK_Cor, fill = id.cluster))
g <- g + geom_violin(trim = FALSE, alpha = 0.7) 
g <- g + geom_boxplot(outlier.shape = NA, width = 0.1,  alpha = 0.7)  +
  theme_classic() +
  ylab("Cytoplasmic log[pERK]") + 
  theme(text = element_text(size = 24),  
        axis.text.x = element_text(angle = 90)) + 
  scale_fill_manual(values = col) + scale_color_manual(values = col)
print(g)

ggsave(filename = "figures/Fig2c_pERK_lineage.png", g,
       width = 450, height = 450, dpi = 300, units = "px", scale = 4)
