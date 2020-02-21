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

theme_set(theme_bw())



# Set working directory to source file


source("../bin/1_Filter_otu_table.R")

############################################# Beta diversity plots and tests 

# Subset: It uses as a base the data set that was created out of the filtered phyloseq object (neg control, relative abundance, no OTUs or Samples with 0 nor 1 abundance) that was later transformed to binary table 


# Use binary_table_OTU_2 (binary table where OTUs that had a presence of 1/ut of all samples, were removed) 

#By site: Texcoco. 
subset.texcoco.binary.beta <- subset_samples(binary_table_OTU2, Site%in%c("mixed", "perturbated", "native"))
subset.texcoco.binary.beta
sample_data(subset.texcoco.binary.beta)

#By site Izta 
subset.izta.binary.beta <- subset_samples(binary_table_OTU2, Site%in%c("Joya", "Cueva", "Base"))
subset.izta.binary.beta
sample_data(subset.izta.binary.beta)

#NMDS Bray Texcoco
bray_nmds = distance(subset.texcoco.binary.beta, method = "bray")
ordination = ordinate(subset.texcoco.binary.beta, method = "NMDS", distance = bray_nmds)
p1 <- plot_ordination(subset.texcoco.binary.beta, ordination, color="Species", shape = "Type", title = "") + theme(aspect.ratio=1)+geom_point(size=3)
print(p1)
p1 + facet_wrap(~Site)

#NMDS Raup-Brick Texcoco

#why using raup-brick: Raup–Crick dissimilarity (method = "raup") is a probabilistic index based on presence/absence data. It is defined as 1 - prob(j), or based on the probability of observing at least j species in shared in compared communities. The current function uses analytic result from hypergeometric distribution (phyper) to find the probabilities. This probability (and the index) is dependent on the number of species missing in both sites, and adding all-zero species to the data or removing missing species from the data will influence the index. The probability (and the index) may be almost zero or almost one for a wide range of parameter values. The index is nonmetric: two communities with no shared species may have a dissimilarity slightly below one, and two identical communities may have dissimilarity slightly above zero. The index uses equal occurrence probabilities for all species, but Raup and Crick originally suggested that sampling probabilities should be proportional to species frequencies (Chase et al. 2011). A simulation approach with unequal species sampling probabilities is implemented in raupcrick function following Chase et al. (2011). The index can be also used for transposed data to give a probabilistic dissimilarity index of species co-occurrence (identical to Veech 2013).
raup_nmds = distance(subset.texcoco.binary.beta, method = "raup")
ordination2 = ordinate(subset.texcoco.binary.beta, method = "NMDS", distance = raup_nmds)
p2 <- plot_ordination(subset.texcoco.binary.beta, ordination2, color="Species", shape = "Type", title = "") + theme(aspect.ratio=1)+geom_point(size=3)
print(p2)
p2 + facet_wrap(~Site)

#NMDS Bray Izta 
bray_nmds2 = distance(subset.izta.binary.beta, method = "bray")
ordination3 = ordinate(subset.izta.binary.beta, method = "NMDS", distance = bray_nmds2)
p3 <- plot_ordination(subset.izta.binary.beta, ordination3, color="Type", title = "") + theme(aspect.ratio=1)+geom_point(size=3)
print(p3)
p3 + facet_wrap(~Site)

#NMDS Raup-Brick Izta
raup_nmds2 = distance(subset.izta.binary.beta, method = "raup")
ordination4 = ordinate(subset.izta.binary.beta, method = "NMDS", distance = raup_nmds2)
p4 <- plot_ordination(subset.izta.binary.beta, ordination4, color="Type", title = "") + theme(aspect.ratio=1)+geom_point(size=3)
print(p4)
p4 + facet_wrap(~Site)

#Community composition Test Texcoco

#Adonis
sampledf <- data.frame(sample_data(subset.texcoco.binary.beta))
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
sampledf2 <- data.frame(sample_data(subset.izta.binary.beta))
adonis2 ( bray_nmds2 ~ Site, data = sampledf2)
adonis2 ( bray_nmds2 ~ Type, data = sampledf2)
adonis2 ( bray_nmds2 ~ Site + Type, data = sampledf2)
adonis2 ( bray_nmds2 ~ Site + (Type/Site), data = sampledf2) #???

adonis2( raup_nmds2 ~ Site, data = sampledf2)
adonis2( raup_nmds2 ~ Type, data = sampledf2)
adonis2( raup_nmds2 ~ Site + (Type/Site), data = sampledf2) #??? 


