#------------------------------------------------------------------------------#
#  This script accompanies the manuscript                                      #
#  Simon et al., (2025)                                                        #
#  Repository available on:                                                    # 
#  https://ctr.uniofcam.dev/ctr-bioinformatics/niakan-lab/simon-et-al-2025     #
#  Please consult READ_ME for more information                                 #
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#                    Control vs FGF Treatments: Analysis                       #
#------------------------------------------------------------------------------#

#..............................................................................#
# Input:    Embryos_clust from FGF_data.R                                      #
#                                                                              #
# Outputs:  Supp. Fig. 1b                                                      #
#           Fig. 1b                                                            #
#           Fig. 1c                                                            #
#           Supp. Fig. 1c                                                      #
#           Supp. Fig. 1d                                                      #
#           Supp. Fig. 1e                                                      #
#           Fig. 1d                                                            #
#           Supp. Fig. 1f                                                      #
#..............................................................................#

# Load library packages, data and set wd ----
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)

setwd("/Users/simonc/Documents/Human_Embryos/Simon_et_al_2025/R")

source("FGF_data.R")


# Colour palettes ----
cols <- c("EPI" = "green", "PrE" = "magenta")
colsTE <- c("EPI" = "green", "PrE" = "magenta",  "TE" = "light blue")
colsICM <- c("ICM" = "grey",  "TE" = "light blue")

# Data wrangling, sample filtering ----

# Count number of cells in each lineage per embryo
Embryos_count <- Embryos_clust %>%
  group_by(Embryo_ID, id.cluster) %>%
  summarize(N = n()) %>%
  ungroup() %>%
  complete(Embryo_ID, id.cluster,
           fill = list(N = 0))

# Re-code treatment labels for each experiment
c25 <- Embryos_count %>% filter(((Embryo_ID == "FGF_1_220318") | 
                                   (Embryo_ID == "FGF_2_220318") ) | 
                                  grepl("250", Embryo_ID)) %>% 
  mutate(Treatment = "250ng/ml FGF")

c5 <- Embryos_count %>% filter(grepl("500", Embryo_ID)) %>% 
  mutate(Treatment = "500ng/ml FGF")

c75 <- Embryos_count %>% filter(grepl("750", Embryo_ID)) %>% 
  mutate(Treatment = "750ng/ml FGF")

cc <- Embryos_count %>% filter(grepl("Control", Embryo_ID)) %>% 
  mutate(Treatment = "Control")

Embryos_count <- rbind(cc, c25, c5, c75)

# Total cell counts per embryo
CellNo <- Embryos_clust %>%
  group_by(Embryo_ID) %>%
  summarize(CellNumber = n())

# Count number of ICM cells per embryo
icmcount <- Embryos_count %>% filter(id.cluster !="TE") %>%
  group_by(Embryo_ID) %>%
  summarize(ICM.No = sum(N)) %>%
  ungroup() %>%
  complete(Embryo_ID,
           fill = list(N = 0))

# Combine counts and calculate percentages in ICM
Embryos_count <- left_join(Embryos_count, icmcount, by = "Embryo_ID")

Embryos_count$perICM <- Embryos_count$N / Embryos_count$ICM.No

Embryos_count <- Embryos_count %>% left_join(CellNo, by = "Embryo_ID")

# Filter out any unfit embryos without an ICM
Embryos_count <- Embryos_count %>% filter(ICM.No != 0)

# Label batches of embryos, experiments before and after revisions
Embryos_count$Batch <- "Batch 1"
Embryos_count <- Embryos_count %>%
  mutate(Batch = ifelse(grepl("2412", Embryo_ID), "Batch 2", Batch))

# Fit embryos with an ICM
withICM <- Embryos_count$Embryo_ID

# Counts of ICM numbers
Embryos_icmcount <- Embryos_clust %>%
  group_by(Embryo_ID, TE_ICM) %>%
  summarize(N = n()) %>%
  ungroup() %>%
  complete(Embryo_ID, TE_ICM,
           fill = list(N = 0))

# Re-code treatment labels for each experiment
ic25 <- Embryos_icmcount %>% 
  filter(((Embryo_ID == "FGF_1_220318") | (Embryo_ID == "FGF_2_220318") ) | 
           grepl("250", Embryo_ID)) %>% mutate(Treatment = "250ng/ml FGF")

ic5 <- Embryos_icmcount %>% 
  filter(grepl("500", Embryo_ID)) %>% mutate(Treatment = "500ng/ml FGF")

ic75 <- Embryos_icmcount %>% filter(grepl("750", Embryo_ID)) %>% 
  mutate(Treatment = "750ng/ml FGF")

icc <- Embryos_icmcount %>% filter(grepl("Control", Embryo_ID)) %>% 
  mutate(Treatment = "Control")

Embryos_icmcount <- rbind(icc, ic25, ic5, ic75)

# Filter out unfit embryos, no ICM
Embryos_icmcount <- Embryos_icmcount %>% filter(Embryo_ID %in% withICM)

# Sample number (n) ----

Embryos_count %>% 
  group_by(Treatment) %>%
  summarise(Unique_Embryo_Count = n_distinct(Embryo_ID))

# # A tibble: 4 × 2
# Treatment    Unique_Embryo_Count
# <chr>                      <int>
# 1 250ng/ml FGF                   8
# 2 500ng/ml FGF                   9
# 3 750ng/ml FGF                   9
# 4 Control                       11

# Supp. Fig 1b: ICM Clustering ----

# Scatter plot of NANOG and GATA4 nuclear fluorescence intensity
# Check hierarchical clustering assignment to epiblast and hypoblast
# Batch 1 and Batch 2 experiments, original submission and revision
# levels changed due to replacement of laser lines on 
# confocal in intervening 2 years

g <- Embryos_clust %>% filter(id.cluster != "TE") %>%
  mutate(Treatment = factor(Treatment, 
                            levels=c("Control", "250ng/ml FGF", 
                                     "500ng/ml FGF", "750ng/ml FGF"))) %>%
  ggplot(aes(y = NANOG_Cor, x = GATA4_Cor))
g <- g + geom_point(aes(color = id.cluster, shape = Treatment), 
                    size = 7, alpha = 0.3)
g <- g + scale_colour_manual(values = cols, 
                             labels = c("Epiblast", "Hypoblast")) + 
  theme_classic() +
  xlab("log[GATA4]") + ylab("log[NANOG]") + 
  theme(text = element_text(size = 32),
        legend.position = "bottom",
        legend.box = "vertical",
        aspect.ratio = 1) + guides(color=guide_legend(title="Lineage")) + 
  facet_wrap(~batch)
print(g)

ggsave(filename = "figures/SuppFig1b_FGF_clustering.png", g,
       width = 4000, height = 2500, dpi = 300, units = "px", scale = 1)



# Fig 1b: Number of hypoblast cells ----

# Boxplot number of hypoblast (PrE) cells in each FGF treatment

g <- Embryos_count %>% filter(id.cluster == "PrE") %>%
  mutate(Treatment = factor(Treatment, 
                            levels=c("Control", "250ng/ml FGF", "500ng/ml FGF", 
                                     "750ng/ml FGF"))) %>%
  ggplot(aes(x = Treatment, y = N))
g <- g + geom_boxplot(aes(fill = id.cluster), outlier.shape = NA, lwd=0.2) 
g <- g + geom_jitter(aes(fill = id.cluster), width = 0, size = 0.8, pch = 21)
g <- g +  stat_compare_means(label = "p.signif", method = "t.test",
                             ref.group = "Control", label.y = 25, size = 3) 
g <- g  + scale_fill_manual(values = cols)  + theme_classic() +
  ylab("Number of hypoblast cells") + theme(text = element_text(size = 8)) +  
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        aspect.ratio = 1)+coord_cartesian(ylim = c(0, 27))
print(g)

ggsave(filename = "figures/Fig1b_FGF_hypoblast.png", g,
       width = 500, height = 650, dpi = 300, units = "px", scale = 1)

# Statistical analysis

# Unpaired t-test comparing mean hypoblast numbers in FGF treatments vs Control

# Control mean = 8

# 250ng/ml FGF hypoblast vs Control 
p25 <- Embryos_count %>% 
  filter(id.cluster == "PrE", 
         Treatment == "Control" | Treatment == "250ng/ml FGF")
t.test(N ~ Treatment, data = p25) # mean = 10  p-value = 0.4883

# 500ng/ml FGF hypoblast vs Control 
p50 <- Embryos_count %>% 
  filter(id.cluster == "PrE", 
         Treatment == "Control" | Treatment == "500ng/ml FGF")
t.test(N ~ Treatment, data = p50) # mean = 13 p-value = 0.1236

# 750ng/ml FGF hypoblast vs Control 
p75 <- Embryos_count %>% 
  filter(id.cluster == "PrE", 
         Treatment == "Control" | Treatment == "750ng/ml FGF")
t.test(N ~ Treatment, data = p75) # mean = 12 p-value = 0.2156

# Fig 1c: Number of epiblast cells ----

# Boxplot number of epiblast cells in each FGF treatment

g <- Embryos_count %>% filter(id.cluster == "EPI") %>%
  mutate(Treatment = factor(Treatment, 
                            levels=c("Control", "250ng/ml FGF", "500ng/ml FGF", 
                                     "750ng/ml FGF"))) %>%
  ggplot(aes(x = Treatment, y = N))
g <- g + geom_boxplot(aes(fill = id.cluster), outlier.shape = NA, lwd=0.2) 
g <- g + geom_jitter(aes(fill = id.cluster), width = 0, size = 0.8, pch = 21)
g <- g +  stat_compare_means(label = "p.signif", method = "t.test",
                             ref.group = "Control", label.y = 20, size = 3) 
g <- g  + scale_fill_manual(values = cols)  + theme_classic() +
  ylab("Number of epiblast cells") + theme(text = element_text(size = 8)) +  
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        aspect.ratio = 1)+coord_cartesian(ylim = c(0, 22))
print(g)

ggsave(filename = "figures/Fig1c_FGF_epiblast.png", g,
       width = 500, height = 650, dpi = 300, units = "px", scale = 1)

# Statistical analysis

# Unpaired t-test comparing mean epiblast numbers in FGF treatments vs Control

# Control mean = 8

# 250ng/ml FGF epiblast vs Control 
e25 <- Embryos_count %>% 
  filter(id.cluster == "EPI", 
         Treatment == "Control" | Treatment == "250ng/ml FGF")
t.test(N ~ Treatment, data = e25) #  mean = 8  p-value = 0.9605

# 500ng/ml FGF epiblast vs Control 
e50 <- Embryos_count %>% 
  filter(id.cluster == "EPI", 
         Treatment == "Control" | Treatment == "500ng/ml FGF")
t.test(N ~ Treatment, data = e50) # mean = 6 p-value =  0.4272

# 750ng/ml FGF epiblast vs Control 
e75 <- Embryos_count %>% 
  filter(id.cluster == "EPI", 
         Treatment == "Control" | Treatment == "750ng/ml FGF")
t.test(N ~ Treatment, data = e75) # mean = 3 p-value = 0.01519



# Supp. Fig. 1c: ICM cell numbers ----

# Boxplot number of ICM cells in each FGF treatment

g <- Embryos_icmcount %>% filter(TE_ICM == "ICM") %>%
  mutate(Treatment = factor(Treatment, 
                            levels=c("Control", "250ng/ml FGF", 
                                     "500ng/ml FGF", "750ng/ml FGF"))) %>%
  ggplot(aes(x = Treatment, y = N, fill = TE_ICM))
g <- g + geom_boxplot(aes(fill = TE_ICM), outlier.shape = NA, lwd=0.2) 
g <- g + geom_jitter(aes(fill = TE_ICM), width = 0, size = 0.8, pch = 21)
g <- g +  stat_compare_means(label = "p.signif", method = "t.test",
                             ref.group = "Control", label.y = 38, size = 3) 
g <- g  + scale_fill_manual(values = colsICM)  + theme_classic() +
  ylab("Number of ICM cells") + theme(text = element_text(size = 8)) +  
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        aspect.ratio = 1)+coord_cartesian(ylim = c(0, 40))
print(g)

ggsave(filename = "figures/SuppFig1c_FGF_ICM.png", g,
       width = 500, height = 650, dpi = 300, units = "px", scale = 1)


# Statistical analysis

# Unpaired t-test comparing mean ICM numbers in FGF treatments vs Control

# Control mean = 15

# 250ng/ml FGF ICM vs Control 
i25 <- Embryos_icmcount %>% 
  filter(TE_ICM == "ICM", Treatment == "Control" | Treatment == "250ng/ml FGF")
t.test(N ~ Treatment, data = i25) # mean = 17 p-value = 0.6582

# 500ng/ml FGF ICM vs Control 
i50 <- Embryos_icmcount %>% 
  filter(TE_ICM == "ICM", Treatment == "Control" | Treatment == "500ng/ml FGF")
t.test(N ~ Treatment, data = i50) # mean = 19 p-value = 0.4971

# 750ng/ml FGF ICM vs Control 
i75 <- Embryos_icmcount %>% filter(TE_ICM == "ICM", Treatment == "Control" | Treatment == "750ng/ml FGF")
t.test(N ~ Treatment, data = i75) # mean = 16 p-value = 0.9482

# Supp. Fig. 1d: TE cell numbers ----

# Boxplot number of TE cells in each FGF treatment

g <- Embryos_count %>% filter(id.cluster == "TE") %>%
  mutate(Treatment = factor(Treatment, 
                            levels=c("Control", "250ng/ml FGF", 
                                     "500ng/ml FGF", "750ng/ml FGF"))) %>%
  ggplot(aes(x = Treatment, y = N))
g <- g + geom_boxplot(aes(fill = id.cluster), outlier.shape = NA, lwd=0.2) 
g <- g + geom_jitter(aes(fill = id.cluster), width = 0, size = 0.8, pch = 21)
g <- g +  stat_compare_means(label = "p.signif", method = "t.test",
                             ref.group = "Control", label.y = 350, size = 3) 
g <- g  + scale_fill_manual(values = colsTE)  + theme_classic() +
  ylab("Number of TE cells") + theme(text = element_text(size = 8)) +  
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        aspect.ratio = 1)+coord_cartesian(ylim = c(0, 400))
print(g)

ggsave(filename = "figures/SuppFig1d_FGF_TE.png", g,
       width = 450, height = 650, dpi = 300, units = "px", scale = 1)



# Supp. Fig. 1e: Total cell numbers ----

# Boxplot number of total cells in each FGF treatment

g <- Embryos_count %>%
  mutate(Treatment = factor(Treatment, levels=c("Control", "250ng/ml FGF", "500ng/ml FGF", "750ng/ml FGF"))) %>%
  ggplot(aes(x = Treatment, y = CellNumber))
g <- g + geom_boxplot(aes(fill = id.cluster), fill="white", outlier.shape = NA, lwd=0.2) 
g <- g + geom_jitter(aes(fill = id.cluster),  fill="white", width = 0, size = 0.8, pch = 21)
g <- g +  stat_compare_means(label = "p.signif", method = "t.test",
                             ref.group = "Control", label.y = 380, size = 3) 
g <- g + theme_classic() +
  ylab("Number of total cells") + theme(text = element_text(size = 8)) +  
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        aspect.ratio = 1)+coord_cartesian(ylim = c(0, 400))
print(g)

ggsave(filename = "figures/SuppFig1e_FGF_cellno.png", g,
       width = 450, height = 650, dpi = 300, units = "px", scale = 1)

# Fig 1d: ICM composition ----

# Stacked barchart of % ICM Cells in each FGF treatment

g <- Embryos_count %>% filter(id.cluster != "TE") %>% 
  mutate(Treatment = factor(Treatment, 
                            levels=c("Control", "250ng/ml FGF", 
                                     "500ng/ml FGF", "750ng/ml FGF"))) %>%
  ggplot(aes(x = Treatment, y = perICM, fill = id.cluster))
g <- g + geom_bar(position="fill", stat="identity")
g <- g + scale_fill_manual(values = cols) + theme_classic() + 
  scale_y_continuous(labels = scales::percent) +
  ylab("% of ICM Cells") + theme(text = element_text(size = 8)) + 
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        aspect.ratio = 0.82)
print(g)

ggsave(filename = "figures/Fig1d_FGF_ICMratio.png", g,
       width = 570, height = 500, dpi = 300, units = "px", scale = 1)

# Re-plot of 1c overlayed with stats (t-test), 
# added to Fig 1c in illustrator manually

g <- Embryos_count %>% filter(id.cluster == "PrE") %>% 
  mutate(Treatment = factor(Treatment, 
                            levels=c("Control", "250ng/ml FGF", 
                                     "500ng/ml FGF", "750ng/ml FGF"))) %>%
  ggplot(aes(x = Treatment, y = perICM, fill = id.cluster))
g <- g + geom_bar(position="fill", stat="identity")
g <- g + scale_fill_manual(values = cols) + theme_classic() + 
  scale_y_continuous(labels = scales::percent) +
  ylab("% of ICM Cells") + theme(text = element_text(size = 8)) + 
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        aspect.ratio = 1)
g <- g +  stat_compare_means(label = "p.signif", method = "t.test",
                             ref.group = "Control") 
print(g)

# Statistical analysis

# Unpaired t-test comparing % hypoblast numbers in FGF treatments vs Control
# Note the p-values are the same when comparing epiblast proportion

# Control 39% PrE

# 250ng/ml FGF hypoblast vs Control 
p25 <- Embryos_count %>% 
  filter(id.cluster == "PrE", 
         Treatment == "Control" | Treatment == "250ng/ml FGF")
t.test(perICM ~ Treatment, data = p25) # p-value = 0.02984 60%

# 500ng/ml FGF hypoblast vs Control 
p50 <- Embryos_count %>% 
  filter(id.cluster == "PrE", 
         Treatment == "Control" | Treatment == "500ng/ml FGF")
t.test(perICM ~ Treatment, data = p50) # p-value = 0.001587 73%

# 750ng/ml FGF hypoblast vs Control 
p75 <- Embryos_count %>% 
  filter(id.cluster == "PrE", 
         Treatment == "Control" | Treatment == "750ng/ml FGF")
t.test(perICM ~ Treatment, data = p75) # p-value = 0.001376 75%

# Supp. Fig. 1f: ICM composition per embryo ----

# Stacked barchart of % ICM Cells for each embryo in each FGF treatment group

g <- Embryos_count %>% filter(id.cluster != "TE") %>% 
  mutate(Treatment = factor(Treatment, 
                            levels=c("Control", "250ng/ml FGF", 
                                     "500ng/ml FGF", "750ng/ml FGF"))) %>%
  ggplot(aes(x = fct_reorder(Embryo_ID, CellNumber), y = perICM, 
             fill = factor(id.cluster, levels = c("EPI", "PrE"))))
g <- g + geom_bar(position="fill", stat="identity")
g <- g + scale_fill_manual(values = cols) + theme_classic() + 
  facet_wrap(~Treatment, drop=TRUE, scales="free_x", nrow = 1) +
  scale_y_continuous(labels = scales::percent) +
  ylab("% of ICM Cells") + theme(legend.position="none",
                                       text = element_text(size = 8),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 aspect.ratio = 1.2)
print(g)

ggsave(filename = "figures/SuppFig1f_FGF_ICMratio_emrbyo.png", g,
       width = 1200, height = 600, dpi = 300, units = "px", scale = 1)
