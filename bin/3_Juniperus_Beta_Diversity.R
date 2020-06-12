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

# Subset by project

# Use binary_table_OTU_2 (binary table where OTUs present in only one sample were removed) 

#Texcoco 
subset.texcoco.binary.beta <- subset_samples(binary_table_OTU2, Project %in% "Texcoco")
subset.texcoco.binary.beta
sample_data(subset.texcoco.binary.beta)

# Remove OTUs that are not present in any of the samples of the subset (OTUs that are only present either in Izta or Texcoco) 
any(taxa_sums(subset.texcoco.binary.beta) == 0)
taxa_sums(subset.texcoco.binary.beta) [1:10]

subset.texcoco.binary.beta <- prune_taxa(taxa_sums(subset.texcoco.binary.beta) > 0, subset.texcoco.binary.beta)

any(taxa_sums(subset.texcoco.binary.beta) == 0)
taxa_sums(subset.texcoco.binary.beta) [1:10]
subset.texcoco.binary.beta


#Izta 
subset.izta.binary.beta <- subset_samples(binary_table_OTU2, Project %in% "Izta")
subset.izta.binary.beta
sample_data(subset.izta.binary.beta)

# Remove OTUs that are not present in any of the samples of the subset (OTUs that are only present either in Izta or Texcoco) 
any(taxa_sums(subset.izta.binary.beta) == 0)
taxa_sums(subset.izta.binary.beta) [1:10]

subset.izta.binary.beta <- prune_taxa(taxa_sums(subset.izta.binary.beta) > 0, subset.izta.binary.beta)

any(taxa_sums(subset.izta.binary.beta) == 0)
taxa_sums(subset.izta.binary.beta) [1:10]
subset.izta.binary.beta


# Data vizualisation

# Reorder Texcoco for facet_wrap
sample_data(subset.texcoco.binary.beta)$Site = factor(sample_data(subset.texcoco.binary.beta)$Site, levels=c("native","mixed","perturbated"))

# NMDS Bray Texcoco
bray_nmds = distance(subset.texcoco.binary.beta, method = "bray")
ordination = ordinate(subset.texcoco.binary.beta, method = "NMDS", distance = bray_nmds)
p1 <- plot_ordination(subset.texcoco.binary.beta, ordination, color="Host", shape = "Type", title = "") + theme(aspect.ratio=1)+geom_point(size=3)
print(p1)
p1 + facet_wrap(~Site)


#NMDS Raup-Brick Texcoco
raup_nmds = distance(subset.texcoco.binary.beta, method = "raup")
ordination2 = ordinate(subset.texcoco.binary.beta, method = "NMDS", distance = raup_nmds)
p2 <- plot_ordination(subset.texcoco.binary.beta, ordination2, color="Host", shape = "Type", title = "") + theme(aspect.ratio=1)+geom_point(size=3)
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




#Community composition Test: Texcoco

#Adonis

sampledf <- data.frame(sample_data(subset.texcoco.binary.beta))
raup_nmds = distance(subset.texcoco.binary.beta, method = "raup")

adonis2( raup_nmds ~ Host, data = sampledf) # significant differences
adonis2( raup_nmds ~ Type, data = sampledf) 
adonis2( raup_nmds ~ Site, data = sampledf) # significant differences 
adonis2( raup_nmds ~ Host + Site + Type, data = sampledf) # significant differences between Hosts and between Sites 

# subset by host plant 
subset.host <- subset_samples(subset.texcoco.binary.beta, Host %in% "Juniperus")
subset.host <- prune_taxa(taxa_sums(subset.host) > 0, subset.host) ## remove OTU with 0 reads
any(taxa_sums(subset.host) == 0)

raup_nmds = distance(subset.host, method = "raup")
sampledf <- data.frame(sample_data(subset.host))

adonis2( raup_nmds ~ Type, data = sampledf)
adonis2( raup_nmds ~ Site, data = sampledf)
adonis2( raup_nmds ~ Site + Type, data = sampledf)
adonis2( raup_nmds ~ Site * Type, data = sampledf)


# subset by site (mixed)
subset.site <- subset_samples(subset.texcoco.binary.beta, Site %in% "mixed")
subset.site <- prune_taxa(taxa_sums(subset.site) > 0, subset.site) ## remove OTU with 0 reads
any(taxa_sums(subset.site) == 0)

raup_nmds = distance(subset.site, method = "raup")
sampledf <- data.frame(sample_data(subset.site))

adonis2( raup_nmds ~ Host, data = sampledf)
adonis2( raup_nmds ~ Type, data = sampledf)
adonis2( raup_nmds ~ Host + Type, data = sampledf)
adonis2( raup_nmds ~ Host * Type, data = sampledf)



adonis2 ( bray_nmds ~ Host + (Type/Host), data = sampledf)
adonis2 ( bray_nmds ~ Host + (Site/Host), data = sampledf) #???? not working because some species are only present in some sites??? 



# Community composition Test Izta

#Adonis
sampledf2 <- data.frame(sample_data(subset.izta.binary.beta))
adonis2 ( bray_nmds2 ~ Site, data = sampledf2) # significant differences 
adonis2 ( bray_nmds2 ~ Type, data = sampledf2) # significant differences
adonis2 ( bray_nmds2 ~ Site + Type, data = sampledf2) #significant differences for both 
adonis2 ( bray_nmds2 ~ Site + (Type/Site), data = sampledf2) # significant differences for both separately but no differences when together

adonis2( raup_nmds2 ~ Site, data = sampledf2) # significant differences 
adonis2( raup_nmds2 ~ Type, data = sampledf2) # significant differences but 0.09 
adonis2( raup_nmds2 ~ Site + (Type/Site), data = sampledf2) # significant differences separately for both, but not together


