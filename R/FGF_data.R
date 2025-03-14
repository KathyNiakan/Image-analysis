#------------------------------------------------------------------------------#
#  This script accompanies the manuscript                                      #
#  Simon et al., (2025)                                                        #
#  Repository available on:                                                    # 
#  https://ctr.uniofcam.dev/ctr-bioinformatics/niakan-lab/simon-et-al-2025     #
#  Please consult READ_ME for more information                                 #
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#                      Control vs FGF Treatments: Data                         #
#------------------------------------------------------------------------------#

#..............................................................................#
# Inputs:   FGF_batch1.csv                                                     #
#           FGF_batch2.csv                                                     #
#           FGF_hclust_batch1.R                                                #
#           FGF_hclust_batch2.R                                                #
# Output:   Embryos_clust                                                      #
#..............................................................................#

# Load library packages and set wd ----
library(dplyr)
library(ggplot2)
library(tidyverse)

setwd("/Users/simonc/Documents/Human_Embryos/Simon_et_al_2025/R")

# Load data ----
Embryos_batch1 <- read_csv("../data/FGF_batch1.csv")
Embryos_batch2 <- read_csv("../data/FGF_batch2.csv")

# Run clustering on ICM cells to separate out Epiblast and Hypoblast ----
source("FGF_hclust_batch1.R")
source("FGF_hclust_batch2.R")

# Merge batches for plotting using FGF_analysis.R
Embryos_clust <- rbind(Embryos_batch1_clust, Embryos_batch2_clust)