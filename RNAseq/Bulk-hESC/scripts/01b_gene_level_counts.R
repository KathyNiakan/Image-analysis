
#module load R/3.5
library(tximport)
library(readr)

# Read in the reference sample data (this file is in folder 'report' on GitHub)
sams <- read_csv("metadata/new_sample_features.csv")

# Subset to new data 
sams[sams$batch == "Simon",]$batch <- "Claire"

# Path to the quantification files
files <- paste("data/tx_quantification", sams$sample_name, "quant.sf", 
               sep = "/")
names(files) <- sams$sample_name

# Load the transcript to gene ID map
load("data/ref_transcriptome/t2g_light.RData")

# Perform the gene-level quantification
all(file.exists(files))
hesc_new <- tximport(files, type = "salmon", tx2gene = t2g)

save(hesc_new, file = "results/hesc_new.RData")
