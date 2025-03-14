#------------------------------------------------------------------------------#
#  This script accompanies the manuscript                                      #
#  Simon et al., (2025)                                                        #
#  Repository available on:                                                    # 
#  https://ctr.uniofcam.dev/ctr-bioinformatics/niakan-lab/simon-et-al-2025     #
#  Please consult READ_ME for more information                                 #
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#                               pERK Human embryos                             #
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
te <- subset(Embryos, TE_ICM == 'TE')
icm <- subset(Embryos, TE_ICM == 'ICM')

# Assign TE cells to a made up cluster 0
te$id.cluster <- "TE"

# Vector of k values to use below
k <- c(4, 3)

# Perform hierarchical clustering on all ICM cells using average linkage
# and cut tree at 4 clusters (k[1])
my.clusters <- hclust(dist(data.frame(icm$SOX2_Cor, 
                                      icm$OTX2_Cor)), 
                      method = 'average')
# Uncomment to see clustering tree
plot(my.clusters)

# Cut tree at k = 4
icm$id.cluster <- cutree(my.clusters, k[1])


## Uncomment below to see outcome
qplot(SOX2_Cor,  OTX2_Cor,
      data = icm, color = id.cluster) + theme_classic() + scale_color_gradient2(low = 'black', mid = 'green',
                                high = 'yellow', midpoint = (k[1]+1)/2)
  
icm$id.cluster[icm$id.cluster==1] <- "EPI"
icm$id.cluster[icm$id.cluster==2] <- "PrE"
icm$id.cluster[icm$id.cluster==3] <- "PrE"
icm$id.cluster[icm$id.cluster==4] <- "PrE"

centers <- icm %>% group_by(id.cluster) %>% summarise(SOX2_Cor = mean(SOX2_Cor), OTX2_Cor = mean(OTX2_Cor))
centers <- centers[-1]
centers <- data.matrix(centers, rownames.force = NA)

Embryos_clust <- rbind(icm, te)
