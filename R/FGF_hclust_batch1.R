#------------------------------------------------------------------------------#
#  This script accompanies the manuscript                                      #
#  Simon et al., (2025)                                                        #
#  Repository available on:                                                    # 
#  https://ctr.uniofcam.dev/ctr-bioinformatics/niakan-lab/simon-et-al-2025     #
#  Please consult READ_ME for more information                                 #
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#                      Control vs FGF Treatments: Batch 1                      #
#------------------------------------------------------------------------------#

################################################################################
#                            Hierarchical clustering                           #
################################################################################
# Modified from hierarchical clustering in Saiz et al (2020) Elife
# Original code available on 
# https://github.com/nestorsaiz/saiz-et-al_2020/blob/master/notebooks/H-clustering.ipynb

# Separate TE and ICM cells
# We will only perform Hierarchical clustering with the icm subset,
# a much smaller dataset, which will be faster, and is the only data  
# we need to classify anyway

te <- subset(Embryos_batch1, TE_ICM == 'TE')
icm <- subset(Embryos_batch1, TE_ICM == 'ICM')

# Assign TE cells to a made up cluster 0
te$id.cluster <- "TE"

# Vector of k values to use below
k <- c(2, 2)

# Perform hierarchical clustering on all ICM cells using average linkage
# and cut tree at 2 clusters (k[2])
my.clusters <- hclust(dist(data.frame(icm$NANOG_Cor, 
                                      icm$GATA4_Cor)), 
                      method = 'average')
# Plot clustering tree
plot(my.clusters)

# Cut tree at k = 2
icm$id.cluster <- cutree(my.clusters, k[2])


# Plot clustering
qplot(NANOG_Cor,  GATA4_Cor,
      data = icm, color = id.cluster) + theme_classic() + scale_color_gradient2(low = 'black', mid = 'green',
                                high = 'yellow', midpoint = (k[2]+1)/2)
  
icm$id.cluster[icm$id.cluster==1] <- "PrE"
icm$id.cluster[icm$id.cluster==2] <- "EPI"

centers <- icm %>% group_by(id.cluster) %>% summarise(NANOG_Cor = mean(NANOG_Cor), GATA4_Cor = mean(GATA4_Cor))
centers <- centers[-1]
centers <- data.matrix(centers, rownames.force = NA)

Embryos_batch1_clust <- rbind(icm, te)
