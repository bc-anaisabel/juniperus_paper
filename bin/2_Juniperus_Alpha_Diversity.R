# Script for Illumina MiSeq2x300 Juniperus projects: Texcoco and Izta
# Part.2 Alpha diversity and direct quantitative observation/comparison of abundances
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

############################################ Alfa diversity 

# Subset:

# By site Texcoco
subset.texcoco.alfa <- subset_samples(juniperus.rel, Site%in%c("mixed", "perturbated", "native"))
subset.texcoco.alfa
sample_data(subset.texcoco.alfa)

# Plotting: 
# fungal species richness per plant species 
plot_richness(subset.texcoco.alfa,x="Species", color = "Site", shape = "Type", measures=c("Observed", "Fisher", "Shannon"))  + geom_point(size=3)

# fungal species richness per site
plot_richness(subset.texcoco.alfa,x="Site", color = "Species", shape = "Type", measures=c("Observed", "Fisher", "Shannon"))  + geom_point(size=3)

# boxplots: 

# fungal species richness per plant species 
plot_richness(subset.texcoco.alfa,x="Species", color = "Site", shape = "Type", measures=c("Observed", "Fisher", "Shannon"))  + geom_boxplot()

# fungal species richness per site
plot_richness(subset.texcoco.alfa,x="Site", color = "Species", shape = "Type", measures=c("Observed", "Fisher", "Shannon"))  + geom_boxplot()

# fungal species richness per site showing only "Observed": 
plot_richness(subset.texcoco.alfa,x="Site", color = "Species", measures=("Observed"))  + geom_boxplot() +facet_wrap(~Type)

# one plot for each diversity measure:
# observed:
a <- plot_richness(subset.texcoco.alfa,x="Site", color = "Species", measures=("Observed"))  
ab <- a + geom_boxplot(data = a$data, aes(x=Site, y=value, color=Species, alpha=0.1)) 
ab +facet_wrap(~Type)

# shannon:
b <- plot_richness(subset.texcoco.alfa, x = "Site", color = "Species", measures=("Shannon")) 
bc <- b + geom_boxplot(data = b$data, aes(x=Site, y=value, color=Species, alpha=0.1))
bc + facet_wrap(~Type)

# fisher:
c <- plot_richness(subset.texcoco.alfa, x = "Site", color = "Species", measures=("Fisher")) 
cc <- c + geom_boxplot(data = c$data, aes(x=Site, y=value, color=Species, alpha=0.1))
cc + facet_wrap(~Type)

# Create a table that gathers diversity measures to use in statistical tests 
texcocodiversity<-estimate_richness(subset.texcoco.alfa, measures=c("Observed", "Fisher", "Shannon"))

data <- cbind(sample_data(subset.texcoco.alfa), texcocodiversity) #combine metadata with alpha diversity

#Site for observed 
subset.texcoco.alfa.anova <- aov(Observed ~ Site, data =data)
summary(subset.texcoco.alfa.anova)
TukeyHSD(subset.texcoco.alfa.anova)
boxplot(Observed~Site, data = data)

#Site for Fisher 
subset.texcoco.alfa.anova <- aov(Fisher ~ Site, data =data)
summary(subset.texcoco.alfa.anova)
TukeyHSD(subset.texcoco.alfa.anova)
boxplot(Fisher~Site, data = data)

#Site for Shannon
subset.texcoco.alfa.anova <- aov(Shannon ~ Site, data =data)
summary(subset.texcoco.alfa.anova)
TukeyHSD(subset.texcoco.alfa.anova)
boxplot(Shannon~Site, data = data)

#Site/Species for observed 
subset.texcoco.alfa.anova2<-aov(Observed~Species+(Species/Site), data = data)
summary(subset.texcoco.alfa.anova2)
TukeyHSD(subset.texcoco.alfa.anova2)
boxplot(Observed~Species+Site, data = data)

#Site/Species for Fisher 
subset.texcoco.alfa.anova2<-aov(Fisher~Species+(Species/Site), data = data)
summary(subset.texcoco.alfa.anova2)
TukeyHSD(subset.texcoco.alfa.anova2)
boxplot(Fisher~Species+Site, data = data)

#Site/Species for Shannon (this one shows significant differences between Plant Species/ Site)
subset.texcoco.alfa.anova2<-aov(Shannon~Species+(Species/Site), data = data)
summary(subset.texcoco.alfa.anova2)
TukeyHSD(subset.texcoco.alfa.anova2)
boxplot(Shannon~Species+Site, data = data)

#Site/Type of sample for observed (this one shows significant differences between Site/Type of sample)
subset.texcoco.alfa.anova <- aov(Observed ~ Site+(Type/Site), data =data)
summary(subset.texcoco.alfa.anova)
TukeyHSD(subset.texcoco.alfa.anova)
boxplot(Observed~Type+Site, data = data)

#Site/Type of sample for Fisher (this one shows significant differences between Site/Type of sample)
subset.texcoco.alfa.anova <- aov(Fisher ~ Site+(Type/Site), data =data)
summary(subset.texcoco.alfa.anova)
TukeyHSD(subset.texcoco.alfa.anova)
boxplot(Fisher~Type+Site, data = data)

#Site/Type of sample for Shannon (this one shows significant differences between Types of sample and between Site/Type of sample)
subset.texcoco.alfa.anova <- aov(Shannon ~ Site+(Type/Site), data =data)
summary(subset.texcoco.alfa.anova)
TukeyHSD(subset.texcoco.alfa.anova)
boxplot(Shannon~Type+Site, data = data)


#? df<- psmelt(subset.texcoco.alfa) #it might be useful to have this data frame 


#Izta: 

#By site Izta 
subset.izta.alfa <- subset_samples(juniperus.rel, Site%in%c("Joya", "Cueva", "Base"))
subset.izta.alfa
sample_data(subset.izta.alfa)

# fungal species richness per type  of sample 
plot_richness(subset.izta.alfa,x="Type", color = "Site", measures=c("Observed", "Fisher", "Shannon"))  + geom_point(size=3)

#fungal species richness per site
plot_richness(subset.izta.alfa,x="Site", color = "Type", measures=c("Observed", "Fisher", "Shannon"))  + geom_point(size=3)


#Do it with boxplots 

# fungal species richness per type of sample
plot_richness(subset.izta.alfa,x="Type", color = "Site", measures=c("Observed", "Fisher", "Shannon")) + geom_boxplot()

# fungal species richness per site
plot_richness(subset.izta.alfa,x="Site", color = "Type", measures=c("Observed", "Fisher", "Shannon"))  + geom_boxplot()

# fungal species richness per site showing only "Observed":  
plot_richness(subset.texcoco.alfa,x="Site", color = "Species", measures=("Observed"))  + geom_boxplot() +facet_wrap(~Type)

# Create a table that gathers diversity measures to use in statistical tests 

iztadiversity<-estimate_richness(subset.izta.alfa, measures=c("Observed", "Fisher", "Shannon"))

data2 <- cbind(sample_data(subset.izta.alfa), iztadiversity) #combine metadata with alpha diversity

#Site for observed 
subset.izta.alfa.anova <- aov(Observed ~ Site, data =data2)
summary(subset.izta.alfa.anova)
TukeyHSD(subset.izta.alfa.anova)
boxplot(Observed~Site, data = data2)

#Site for Fisher 
subset.izta.alfa.anova <- aov(Fisher ~ Site, data =data2)
summary(subset.izta.alfa.anova)
TukeyHSD(subset.izta.alfa.anova)
boxplot(Fisher~Site, data = data2)

#Site for Shannon
subset.izta.alfa.anova <- aov(Shannon ~ Site, data =data2)
summary(subset.izta.alfa.anova)
TukeyHSD(subset.izta.alfa.anova)
boxplot(Shannon~Site, data = data2)

#Site/Type of sample for observed (this one shows significant differences between types of sample)
subset.izta.alfa.anova <- aov(Observed ~ Site+(Type/Site), data =data2)
summary(subset.izta.alfa.anova)
TukeyHSD(subset.izta.alfa.anova)
boxplot(Observed~Type+Site, data = data2)

#Site/Type of sample for Fisher (this one shows significant differences between types of sample)
subset.izta.alfa.anova <- aov(Fisher ~ Site+(Type/Site), data =data2)
summary(subset.izta.alfa.anova)
TukeyHSD(subset.izta.alfa.anova)
boxplot(Fisher~Type+Site, data = data2)

#Site/Type of sample for Shannon (this one shows significant differences between types of sample)
subset.izta.alfa.anova <- aov(Shannon ~ Site+(Type/Site), data =data2)
summary(subset.izta.alfa.anova)
TukeyHSD(subset.izta.alfa.anova)
boxplot(Shannon~Type+Site, data = data2)

# ? df2<- psmelt(subset.izta.alfa) #it might be useful to have this data frame 



################################## Exploration for direct quantitative observation/comparison of abundances

a<-taxa_sums(juniperus.rel)[1:100]
print(a)
sample_sums(subset.texcoco.alfa) [1:10] #Should we remove OTUs that are not present in any of the samples of the subset? (OTUs that are only present either in Izta or Texcoco) if we do, do it in alpha diversity and beta diversity analysis? 
taxa_sums(subset.texcoco.alfa) [1:10]

# get_taxa is returning all the OTU abundances from one sample, while get_sample is returning the abundances from all samples for one OTU.
get_taxa(subset.texcoco.alfa, sample_names(subset.texcoco.alfa)[5])[1:10]
get_sample(subset.texcoco.alfa, taxa_names(subset.texcoco.alfa)[5])[1:10]

# melt to long format (for ggploting) 
# prune out phyla below 1% in each sample
# selecting the taxa at the level: Phylum

mdata_phylum <- subset.texcoco.alfa %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

# checking the dataframe that we now created
head(mdata_phylum)

#Now plot by Relative Abundance of Phylum by Site, Type of sample and Plant species 
ggplot(mdata_phylum, aes(x = Site, y = Abundance, fill = Phylum)) + 
  #facet_grid(time~.) +
  geom_bar(stat = "identity")  +
  # Remove x axis title, and rotate sample lables
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # additional stuff
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("Relative Abundance (Phyla > 1%)\n") +
  ggtitle("Phylum Relative Abundance") + facet_grid(Type ~ sample_Species)



# Top 50 OTUs
Top50OTUs <- names(sort(taxa_sums(subset.texcoco.alfa), TRUE)[1:50])
ent50 <- prune_taxa(Top50OTUs, subset.texcoco.alfa)
plot_bar(ent50, "Site", fill = "Family")
plot_bar(ent50, "sample_Species", fill = "Family", facet_grid = "Site")

# only Ascomycota, Basidiomycota and Glomeromycota
subset.phylum <- subset_taxa(subset.texcoco.alfa, Phylum %in% c("p__Ascomycota", "p__Basidiomycota", "p__Glomeromycota"))
plot_bar(subset.phylum, "Site", fill = "Phylum", facet_grid = Project ~ sample_Species)

# only Glomeromycota 
subset.phylum.glomero <- subset_taxa(subset.texcoco.alfa, Phylum =="p__Glomeromycota")
plot_bar(subset.phylum.glomero, "Site", fill = "Phylum", facet_grid = Project ~ sample_Species)

# Phylum
plot_bar(subset.texcoco.alfa , "Site", fill = "Phylum") + facet_wrap(sample_Species ~ Type) 

# Family
plot_bar(ent50, fill = "Family")
plot_bar(ent50, fill = "Family") + facet_wrap(~Site, scales="free_x", nrow=1)
plot_bar(ent50, fill = "Family") + facet_wrap(~sample_Species+ Site + Type, scales="free_x", nrow=1)




