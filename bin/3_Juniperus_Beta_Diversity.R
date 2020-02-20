# Script for Illumina MiSeq2x300 Juniperus projects: Texcoco and Izta
# Part. 3. Beta diversity 
# February, 2020 

library("plyr"); packageVersion("plyr")
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("vegan"); packageVersion("vegan")
library("RColorBrewer"); packageVersion("RColorBrewer")
library("plotly"); packageVersion("plotly")
library("htmltools"); packageVersion("htmltools")
library("DT"); packageVersion("DT")
library(ggplot2)
library(dplyr)
library(tibble)


setwd("~/Documents/Ana+Camille/MiSeq2x300_Dec2019/Juniperus/data/taxonomy")

source("~/Documents/Ana+Camille/MiSeq2x300_Dec2019/Juniperus/R/Script/1_Juniperus_Filters.R")

############################################# Beta diversity plots and tests 

#subset:

#By site Texcoco
subset.texcoco <- subset_samples(subset.1, Site%in%c("mixed", "perturbated", "native"))
subset.texcoco
sample_data(subset.texcoco)

#By site Izta 
subset.izta <- subset_samples(subset.1, Site%in%c("Joya", "Cueva", "Base"))
subset.izta
sample_data(subset.izta)

#NMDS Bray Texcoco
bray_nmds = distance(subset.texcoco, method = "bray")
ordination = ordinate(subset.texcoco, method = "NMDS", distance = bray_nmds)
p1 <- plot_ordination(subset.texcoco, ordination, color="Species", shape = "Type", title = "") + theme(aspect.ratio=1)+geom_point(size=3)
print(p1)
p1 + facet_wrap(~Site)

#NMDS Raup-Brick Texcoco

#why using raup-brick: Raup–Crick dissimilarity (method = "raup") is a probabilistic index based on presence/absence data. It is defined as 1 - prob(j), or based on the probability of observing at least j species in shared in compared communities. The current function uses analytic result from hypergeometric distribution (phyper) to find the probabilities. This probability (and the index) is dependent on the number of species missing in both sites, and adding all-zero species to the data or removing missing species from the data will influence the index. The probability (and the index) may be almost zero or almost one for a wide range of parameter values. The index is nonmetric: two communities with no shared species may have a dissimilarity slightly below one, and two identical communities may have dissimilarity slightly above zero. The index uses equal occurrence probabilities for all species, but Raup and Crick originally suggested that sampling probabilities should be proportional to species frequencies (Chase et al. 2011). A simulation approach with unequal species sampling probabilities is implemented in raupcrick function following Chase et al. (2011). The index can be also used for transposed data to give a probabilistic dissimilarity index of species co-occurrence (identical to Veech 2013).
raup_nmds = distance(subset.texcoco, method = "raup")
ordination2 = ordinate(subset.texcoco, method = "NMDS", distance = raup_nmds)
p2 <- plot_ordination(subset.texcoco, ordination2, color="Species", shape = "Type", title = "") + theme(aspect.ratio=1)+geom_point(size=3)
print(p2)
p2 + facet_wrap(~Site)

#NMDS Bray Izta 
bray_nmds2 = distance(subset.izta, method = "bray")
ordination3 = ordinate(subset.izta, method = "NMDS", distance = bray_nmds2)
p3 <- plot_ordination(subset.izta, ordination3, color="Type", title = "") + theme(aspect.ratio=1)+geom_point(size=3)
print(p3)
p3 + facet_wrap(~Site)

#NMDS Raup-Brick Izta
raup_nmds2 = distance(subset.izta, method = "raup")
ordination4 = ordinate(subset.izta, method = "NMDS", distance = raup_nmds2)
p4 <- plot_ordination(subset.izta, ordination4, color="Type", title = "") + theme(aspect.ratio=1)+geom_point(size=3)
print(p4)
p4 + facet_wrap(~Site)

#Community composition Test Texcoco

#Adonis
sampledf <- data.frame(sample_data(subset.texcoco))
adonis2 ( bray_nmds ~ Species + Type, data = sampledf)
adonis2 ( bray_nmds ~ Species + Site, data = sampledf)
adonis2 ( bray_nmds ~ Species + (Type/Species), data = sampledf)
adonis2 ( bray_nmds ~ Species + (Site/Species), data = sampledf) #???? not working because some species are only present in some sites??? 

adonis2( raup_nmds ~ Species + Type, data = sampledf)
adonis2( raup_nmds ~ Species + (Type/Species), data = sampledf)
adonis2( raup_nmds ~ Species + Site, data = sampledf)
adonis2( raup_nmds ~ Species + (Site/Species), data = sampledf) #???? not working because some species are only present in some sites??? 

# Community composition Test Izta

#Adonis
sampledf2 <- data.frame(sample_data(subset.izta))
adonis2 ( bray_nmds2 ~ Site, data = sampledf2)
adonis2 ( bray_nmds2 ~ Type, data = sampledf2)
adonis2 ( bray_nmds2 ~ Site + Type, data = sampledf2)
adonis2 ( bray_nmds2 ~ Site + (Type/Site), data = sampledf2) #???

adonis2( raup_nmds2 ~ Site, data = sampledf2)
adonis2( raup_nmds2 ~ Type, data = sampledf2)
adonis2( raup_nmds2 ~ Site + (Type/Site), data = sampledf2) #??? 


