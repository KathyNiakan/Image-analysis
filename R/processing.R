#  This script accompanies the manuscript                                      #
#  Simon et al., (2025)                                                        #
#  Repository available on:                                                    # 
#  https://ctr.uniofcam.dev/ctr-bioinformatics/niakan-lab/simon-et-al-2025     #
#  Please consult READ_ME for more information                                 #
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#    Processing CellProfiler data: Example 220225 Control 1 and 2 (ERKi)       #
#------------------------------------------------------------------------------#

#..............................................................................#
# Inputs:   Stardist_testObjects_edge_size_filtered.csv                        #
#           Segmentation.csv                                                   #
#           ICM.csv                                                            #  
#           ebcor.R                                                            #
# Output:   Embryos (data frame)                                               #
#..............................................................................#

# Load library packages, set wd and load scripts ----
library(dplyr)
library(ggplot2)
library(tidyverse)

setwd("/Users/simonc/Documents/Human_Embryos/Simon_et_al_2025/R")
source("ebcor.R")

# Load and data wrangling of CellProfiler and manual corrections ----

# Specify the Embryo IDs and therefore parent folders
emb_no <- c("Control1", "Control2")

# Specify path where CellProfiler output is
parent_path <- c("/Users/simonc/Documents/Human_Embryos/Simon_et_al_2025/data/embryo_example")
setwd(parent_path)

# Data wrangling CellProfiler output
# Specify the variables to keep for nuclear objects
NucVars <- c("ImageNumber", "ObjectNumber", "AreaShape_Area", 
             "Intensity_MeanIntensity_C0", "Intensity_MeanIntensity_C1", 
             "Intensity_MeanIntensity_C2", "Intensity_MeanIntensity_C3",
             "Intensity_MeanIntensity_C4", "Location_Center_X", 
             "Location_Center_Y", "Location_Center_Z", "TrackObjects_Label_1",
             "TrackObjects_Lifetime_1")

# # For pERK staining to measure cytoplasmic values uncomment the following lines
# CytoVars <- c("Intensity_MeanIntensity_C2", "IO")
# CytoRename <- c("CytoIntensity_MeanIntensity_C2", "IO")

Embryos <- NULL
  for (i in emb_no) {

Nuc <- read.csv(paste0(parent_path, "/CellProfiler_", i, 
                       "/Stardist_testObjects_edge_size_filtered.csv"))
Nuc$Location_Center_Z <- Nuc$ImageNumber
NucSub <- Nuc[NucVars]

# # For pERK staining to measure cytoplasmic values uncomment the following lines
# Cyto <- read.csv(paste0(parent_path, "/CellProfiler_", i, "/Stardist_testCytoring.csv"))
# Cyto$Location_Center_Z <- Cyto$ImageNumber
# Cyto$IO <- paste0(Cyto$ImageNumber, sep = "_", Cyto$ObjectNumber)
# CytoSub <- Cyto[CytoVars]
# names(CytoSub) <- CytoRename
# All <- left_join(NucSub, CytoSub, by = "IO")


# Rename nuclear intensities for each channel, average across all Z slices,
# get positional data
TrackedCells <- NucSub %>% filter(!is.na(TrackObjects_Label_1)) %>% # For pERK replace NucSub with All
  group_by(TrackObjects_Label_1) %>% 
  summarise(DAPI = mean(Intensity_MeanIntensity_C0)*100, 
            NANOG = mean(Intensity_MeanIntensity_C1)*100, # For pERK exp this is Sox2
            GATA4 = mean(Intensity_MeanIntensity_C2)*100, # For pERK exp this is pERK (nuclear)
            GATA3 = mean(Intensity_MeanIntensity_C3)*100, # For pERK exp this is BF
            BF = mean(Intensity_MeanIntensity_C4)*100, # For pERK exp this is Otx2
           # Cyto_pERK = mean(CytoIntensity_MeanIntensity_C2)*100, # For pERK cytoplasmic uncomment
            X = mean(Location_Center_X),
            Y = mean(Location_Center_Y),
            Z = mean(Location_Center_Z),
            Lifetime = max(TrackObjects_Label_1),
            Area = sum(AreaShape_Area)) 

# Add experiment date to Embryo ID
TrackedCells$Embryo_ID <- paste0(i,"_220205")

# Fetch manual label on which are ICM vs TE cells and add label
ICM_label <- read.csv(paste0(parent_path, "/CellProfiler_", i, "/ICM.csv"))
ICM_ls <- ICM_label$TrackObjects_Label_1

ICM <- TrackedCells %>% group_by(TrackObjects_Label_1) %>% 
  filter(TrackObjects_Label_1 %in% ICM_ls) %>% mutate(TE_ICM = "ICM")
TE <- TrackedCells %>% group_by(TrackObjects_Label_1) %>% 
  filter(!TrackObjects_Label_1 %in% ICM_ls) %>% mutate(TE_ICM = "TE")

# Fetch data with manual corrections for under segmentation
f <- paste0(parent_path, "/CellProfiler_", i, "/Segmentation.csv")
if (file.exists(f)) {
corrected <- read.csv(f)
under <- read.csv(paste0(parent_path, "/CellProfiler_", i, "/ICM.csv"))
under <- under %>% filter(Segmentation == "UNDER")
under <- under$TrackObjects_Label_1
ICM <- ICM %>% filter(!(TrackObjects_Label_1 %in% under))
Embi <- rbind(ICM, TE, corrected)
} else {
Embi <- rbind(ICM, TE) }


# Embi <- rbind(ICM, TE)
Embryos <- rbind(Embryos, Embi)

  }

# Apply Z-correction for nuclear fluorescence intensities ----

# Note: requires n >= 2 embryos to run ebcor script

# Define the possible channels to go through
channels <- c('DAPI', 'NANOG', 'GATA4', 'GATA3', 'BF') # For pERK change channel names and include Cyto_pERK
# Create an empty matrix to hold the corrected values
# with the length of the number of cells in the dataset
ebLogCor <- matrix(0, nrow = length(Embryos$Embryo_ID), 
                   # and as many columns as channels
                   ncol = length(channels),
                   dimnames = list(c(), channels))
# For each channel in 'channels', perform EB correction 
# and store in the corresponding column in 'ebLogCor'
for (c in channels) {
  ebLogCor[, c] <- ebcor(Embryos, c)
}

# Convert ebLogCor to data frame and rename columns to MARKER_Cor
ebLogCor <- data.frame(ebLogCor)
ebLogCor <- rename(ebLogCor, DAPI_Cor = DAPI, 
                   NANOG_Cor = NANOG,
                   GATA4_Cor = GATA4, 
                   GATA3_Cor = GATA3,
                   BF_Cor = BF) # For pERK change channel names and include Cyto_pERK

# Incorporate ebLogCor values into main table
Embryos <- cbind(Embryos, ebLogCor)
rm(ebLogCor)

# Add treatment information and re-bind data
# Uncomment if multiple treatments types
Control <- Embryos %>% filter(grepl("Control", Embryo_ID)) %>% 
  mutate(Treatment = "Control")
# Ulix <- Embryos %>% filter(grepl("Ulix", Embryo_ID)) %>% 
  mutate(Treatment = "Ulix")
Embryos <- Control
# Embryos <- rbind(Control, Ulix)

# Processed data can be saved and analysed using downstream scripts
