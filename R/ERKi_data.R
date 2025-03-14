#------------------------------------------------------------------------------#
#  This script accompanies the manuscript                                      #
#  Simon et al., (2025)                                                        #
#  Repository available on:                                                    # 
#  https://ctr.uniofcam.dev/ctr-bioinformatics/niakan-lab/simon-et-al-2025     #
#  Please consult READ_ME for more information                                 #
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#                      Control vs ERKi Treatment: Data                         #
#------------------------------------------------------------------------------#

#..............................................................................#
# Inputs:   ERKi.csv                                                           #
#           ERKi_hclust.R                                                      #
# Output:   Embryos_clust                                                      #
#..............................................................................#

# Load library packages and set wd ----
library(dplyr)
library(ggplot2)
library(tidyverse)

setwd("/Users/simonc/Documents/Human_Embryos/Simon_et_al_2025/R")

# Load data ----
Embryos <- read_csv("../data/ERKi.csv")

# Run clustering on ICM cells to separate out Epiblast and Hypoblast ----
source("ERKi_hclust.R")