#------------------------------------------------------------------------------#
#  This script accompanies the manuscript                                      #
#  Simon et al., (2025)                                                        #
#  Repository available on:                                                    # 
#  https://ctr.uniofcam.dev/ctr-bioinformatics/niakan-lab/simon-et-al-2025     #
#  Please consult READ_ME for more information                                 #
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#                    Control vs ERKi Treatments: Analysis                      #
#------------------------------------------------------------------------------#

#..............................................................................#
# Input:    Embryos_clust from ERKi_data.R                                     #
#                                                                              #
# Outputs:  Supp. Fig. 2e                                                      #
#           Fig. 2e                                                            #
#           Fig. 2f                                                            #
#           Supp. Fig. 2f                                                      #
#           Fig. 2g                                                            #
#           Supp. Fig. 2g                                                      #
#           Supp. Fig. 2j                                                      #
#           Supp. Fig. 2k                                                      #
#..............................................................................#

# Load library packages, data and set wd ----
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)

setwd("/Users/simonc/Documents/Human_Embryos/Simon_et_al_2025/R")

source("ERKi_data.R")

# Colour palettes and plot labels ----

# Palette
cols <- c("EPI" = "green", "PrE" = "magenta")
colsTE <- c("EPI" = "green", "PrE" = "magenta", "TE" = "light blue")
colsICM <- c("ICM" = "grey",  "TE" = "light blue")

# Labels
treat.labs <- c("Control", "ERKi")
names(treat.labs) <- c("Control", "Ulix")

# Data wrangling, sample filtering ----

# Count number of cells in each lineage per embryo
Embryos_count <- Embryos_clust %>%
  group_by(Embryo_ID, id.cluster) %>%
  summarize(N = n()) %>%
  ungroup() %>%
  complete(Embryo_ID, id.cluster,
           fill = list(N = 0))

# Re-code treatment labels for each experiment
cu <- Embryos_count %>% filter(grepl("Ulix", Embryo_ID)) %>% 
  mutate(Treatment = "Ulix")
cc <- Embryos_count %>% filter(grepl("Control", Embryo_ID)) %>% 
  mutate(Treatment = "Control")
Embryos_count <- rbind(cc, cu)

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

# Counts of ICM numbers
Embryos_icmcount <- Embryos_clust %>%
  group_by(Embryo_ID, TE_ICM) %>%
  summarize(N = n()) %>%
  ungroup() %>%
  complete(Embryo_ID, TE_ICM,
           fill = list(N = 0))

# Re-code treatment labels for each experiment
icu <- Embryos_icmcount %>% filter(grepl("Ulix", Embryo_ID)) %>% 
  mutate(Treatment = "Ulix")

icc <- Embryos_icmcount %>% filter(grepl("Control", Embryo_ID)) %>% 
  mutate(Treatment = "Control")

Embryos_icmcount <- rbind(icc, icu)

# Sample number (n) ----

Embryos_count %>% 
  group_by(Treatment) %>%
  summarise(Unique_Embryo_Count = n_distinct(Embryo_ID))


# Supp. Fig 2e: ICM Clustering ----

# Scatter plot of NANOG and GATA4 nuclear fluorescence intensity
# Check hierarchical clustering assignment to epiblast and hypoblast


g <- Embryos_clust %>% filter(id.cluster != "TE") %>%
  ggplot(aes(y = NANOG_Cor, x = GATA4_Cor))
g <- g + geom_point(aes(color = id.cluster, shape = Treatment), 
                    size = 6, alpha = 0.3)
g <- g + scale_colour_manual(values = cols, 
                             labels = c("Epiblast", "Hypoblast")) + 
  scale_shape_discrete(labels = c("Control", "ERKi"))+ theme_classic() +
  xlab("log[GATA4]") + ylab("log[NANOG]") + 
  theme(text = element_text(size = 34),
        legend.position = "right",
        legend.box = "vertical",
        aspect.ratio = 1) + guides(color=guide_legend(title="Lineage"))
print(g)

ggsave(filename = "figures/SuppFig2_ERKi_clust.png", g,
       width = 2800, height = 2500, dpi = 300, units = "px", scale = 1)


# Fig. 2e: Number of hypoblast cells ----

# Boxplot number of hypoblast (PrE) cells in each FGF treatment

g <- Embryos_count %>% filter(id.cluster == "PrE") %>%
  ggplot(aes(x = Treatment, y = N, fill = id.cluster))
g <- g + geom_boxplot(aes(fill = id.cluster), outlier.shape = NA, lwd=0.2) 
g <- g + geom_jitter(aes(fill = id.cluster), width = 0, size = 0.8, pch = 21)
g <- g +  stat_compare_means(label = "p.signif", method = "t.test",
                             ref.group = "Control", label.y = 28, size = 3) 
g <- g  + scale_fill_manual(values = cols)  + theme_classic() +
  ylab("Number of hypoblast cells") + theme(text = element_text(size = 8)) +  
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        aspect.ratio = 2)+coord_cartesian(ylim = c(0, 30)) +
  scale_x_discrete(labels=c("Ulix" = "ERKi"))
print(g)

ggsave(filename = "figures/Fig2e_ERKi_hypoblast.png", g,
       width = 350, height = 650, dpi = 300, units = "px", scale = 1)


# Statistical analysis

# Unpaired t-test comparing mean hypoblast numbers in ERKi treatment vs Control

# Control mean = 13

# ERKi hypoblast vs Control
PrEc <- Embryos_count %>% filter(id.cluster == "PrE")
t.test(N ~ Treatment, data = PrEc) # mean = 2 p-value = 0.002905


# Fig. 2f: Number of epiblast cells ----

# Boxplot number of epiblast cells in ERKi treatment

g <- Embryos_count %>% filter(id.cluster == "EPI") %>%
  ggplot(aes(x = Treatment, y = N, fill = id.cluster))
g <- g + geom_boxplot(aes(fill = id.cluster), outlier.shape = NA, lwd=0.2) 
g <- g + geom_jitter(aes(fill = id.cluster), width = 0, size = 0.8, pch = 21)
g <- g +  stat_compare_means(label = "p.signif", method = "t.test",
                             ref.group = "Control", label.y = 33, size = 3) 
g <- g  + scale_fill_manual(values = cols)  + theme_classic() +
  ylab("Number of epiblast cells") + theme(text = element_text(size = 8)) +  
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        aspect.ratio = 2)+coord_cartesian(ylim = c(0, 35)) +
  scale_x_discrete(labels=c("Ulix" = "ERKi"))
print(g)

ggsave(filename = "figures/Fig2e_ERKi_epiblast.png", g,
       width = 350, height = 650, dpi = 300, units = "px", scale = 1)


# Statistical analysis

# Unpaired t-test comparing mean epiblast numbers in ERKi treatment vs Control

# Control mean = 12

# ERKi epiblast vs Control
EPIc <- Embryos_count %>% filter(id.cluster == "EPI")

t.test(N ~ Treatment, data = EPIc) # mean = 15 p-value = 0.3996




# Supp. Fig. 2f: ICM cell numbers ----

# Boxplot number of ICM cells in ERKi treatment

g <- Embryos_icmcount %>% filter(TE_ICM == "ICM") %>%
  ggplot(aes(x = Treatment, y = N, fill = TE_ICM))
g <- g + geom_boxplot(aes(fill = TE_ICM), outlier.shape = NA, lwd=0.2) 
g <- g + geom_jitter(aes(fill = TE_ICM), width = 0, size = 0.8, pch = 21)
g <- g +  stat_compare_means(label = "p.signif", method = "t.test",
                             ref.group = "Control", label.y = 48, size = 3) 
g <- g  + scale_fill_manual(values = colsICM)  + theme_classic() +
  ylab("Number of ICM cells") + theme(text = element_text(size = 8)) +  
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        aspect.ratio = 1.5)+coord_cartesian(ylim = c(0, 50)) +
  scale_x_discrete(labels=c("Ulix" = "ERKi"))
print(g)

ggsave(filename = "figures/SuppFig2e_ERKi_ICM.png", g,
       width = 350, height = 600, dpi = 300, units = "px", scale = 1)

# Statistical analysis

# Unpaired t-test comparing mean ICM numbers in ERKi treatment vs Control

# Control mean = 25

# ERKi ICM vs Control
i <- Embryos_icmcount %>% filter(TE_ICM == "ICM")
t.test(N ~ Treatment, data = i) #  mean = 17  p-value = 0.129


# Fig. 2g: ICM composition ----

# Stacked barchart of % ICM Cells in ERKi treatment

g <- Embryos_count %>% filter(id.cluster != "TE") %>%
  ggplot(aes(x = Treatment, y = perICM, fill = id.cluster))
g <- g + geom_bar(position="fill", stat="identity")
g <- g + scale_fill_manual(values = cols) + theme_classic() + 
  scale_y_continuous(labels = scales::percent) +
  ylab("% of ICM Cells") + theme(text = element_text(size = 8)) + 
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        aspect.ratio = 2) +
  scale_x_discrete(labels=c("Ulix" = "ERKi"))

print(g)

ggsave(filename = "figures/Fig2e_ERKi_ICMratio.png", g,
       width = 350, height = 500, dpi = 300, units = "px", scale = 1)

# Re-plot of 2g overlayed with stats (t-test), 
# added to Fig 2g in illustrator manually

g <- Embryos_count %>% filter(id.cluster == "PrE") %>%
  ggplot(aes(x = Treatment, y = perICM, fill = id.cluster))
g <- g + geom_bar(position="fill", stat="identity")
g <- g + scale_fill_manual(values = cols) + theme_classic() + 
  scale_y_continuous(labels = scales::percent) +
  ylab("% of ICM Cells") + theme(text = element_text(size = 8)) + 
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        aspect.ratio = 2) +
  scale_x_discrete(labels=c("Ulix" = "ERKi"))
g <- g +  stat_compare_means(label = "p.signif", method = "t.test",
                             ref.group = "Control") 

print(g)


# Statistical analysis

# Unpaired t-test comparing % hypoblast numbers in FGF treatments vs Control
# Note the p-values are the same when comparing epiblast proportion

# Control 52% PrE

# ERKi hypoblast vs Control 
p <- Embryos_count %>% filter(id.cluster == "PrE")
t.test(perICM ~ Treatment, data = p) # mean = 9%  p-value = 0.0002288


# Supp. Fig. 2g: ICM composition per embryo ----

# Stacked barchart of % ICM Cells for each embryo in each FGF treatment group

g <- Embryos_count %>% filter(id.cluster != "TE") %>% 
  ggplot(aes(x = fct_reorder(Embryo_ID, CellNumber), y = N, 
             fill = factor(id.cluster, levels = c("EPI", "PrE"))))
g <- g + geom_bar(position="fill", stat="identity")
g <- g + scale_fill_manual(values = cols) + theme_classic() + 
  facet_wrap(~Treatment, drop=TRUE, scales="free", nrow = 1, 
             labeller = labeller(Treatment = treat.labs)) +
  scale_y_continuous(labels = scales::percent) +
  ylab("% of ICM Cells") + theme(legend.position="none",
                                 text = element_text(size = 8),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 aspect.ratio = 1.5)
print(g)

ggsave(filename = "figures/SuppFig2g_ERKi_ICMratio_embryo.png", g,
       width = 700, height = 500, dpi = 300, units = "px", scale = 1)


# Supp. Fig. 2j: TE cell numbers ----

# Boxplot number of TE cells in ERKi treatment

g <- Embryos_icmcount %>% filter(TE_ICM == "TE") %>%
  ggplot(aes(x = Treatment, y = N, fill = TE_ICM))
g <- g + geom_boxplot(aes(fill = TE_ICM), outlier.shape = NA, lwd=0.2) 
g <- g + geom_jitter(aes(fill = TE_ICM), width = 0, size = 0.8, pch = 21)
g <- g +  stat_compare_means(label = "p.signif", method = "t.test",
                             ref.group = "Control", label.y = 408, size = 3) 
g <- g  + scale_fill_manual(values = colsICM)  + theme_classic() +
  ylab("Number of TE cells") + theme(text = element_text(size = 8)) +  
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        aspect.ratio = 2)+coord_cartesian(ylim = c(0, 430)) +
  scale_x_discrete(labels=c("Ulix" = "ERKi"))
print(g)

ggsave(filename = "figures/SuppFig2j_ERKi_TE.png", g,
       width = 300, height = 600, dpi = 300, units = "px", scale = 1)


# Statistical analysis

# Unpaired t-test comparing mean TE numbers in ERKi treatment vs Control

# Control mean = 272

# ERKi ICM vs Control
t <- Embryos_icmcount %>% filter(TE_ICM == "TE")
t.test(N ~ Treatment, data = t) # mean = 158,  p-value = 0.03881



# Supp. Fig. 2k: Total cell numbers ----

# Boxplot number of total cells in ERKi treatment

g <- Embryos_count %>%
  ggplot(aes(x = Treatment, y = CellNumber))
g <- g + geom_boxplot(fill = "white", outlier.shape = NA, lwd=0.2) 
g <- g + geom_jitter(fill = "white", width = 0, size = 0.8, pch = 21)
g <- g +  stat_compare_means(label = "p.signif", method = "t.test",
                             ref.group = "Control", label.y = 448, size = 3) 
g <- g + theme_classic() +
  ylab("Number of total cells") + theme(text = element_text(size = 8)) +  
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        aspect.ratio = 2)+coord_cartesian(ylim = c(0, 470)) +
  scale_x_discrete(labels=c("Ulix" = "ERKi"))
print(g)

ggsave(filename = "figures/SuppFig2k_ERKi_cellno.png", g,
       width = 300, height = 600, dpi = 300, units = "px", scale = 1)


# Statistical analysis

# Unpaired t-test comparing mean TE numbers in ERKi treatment vs Control

# Control mean = 297
c <- Embryos_count %>% group_by(Embryo_ID, CellNumber, Treatment) %>% summarise()
t.test(CellNumber ~ Treatment, data = c) # mean = 175 p-value = 0.03881