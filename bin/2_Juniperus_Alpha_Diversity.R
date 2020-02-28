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


# Set working directory to source file


source("../bin/1_Filter_otu_table.R")



################################## Exploration for direct quantitative observation/comparison of abundances

######################## Subset: Texcoco 

# Subset data for Texcoco using OTU table with relative abundance 
subset.texcoco.alfa <- subset_samples(phyloseq.rel, Project %in% "Texcoco")
subset.texcoco.alfa

sample_data(subset.texcoco.alfa)
sample_sums(subset.texcoco.alfa) [1:10] 

# Remove OTUs that are not present in any of the samples of the subset (OTUs that are only present either in Izta or Texcoco) 
taxa_sums(subset.texcoco.alfa) [1:10]

subset.texcoco.alfa <- prune_taxa(taxa_sums(subset.texcoco.alfa) > 0, subset.texcoco.alfa)
any(taxa_sums(subset.texcoco.alfa) == 0)

taxa_sums(subset.texcoco.alfa) [1:10]
subset.texcoco.alfa

# Subset data for Texcoco using OTU binary table 

subset.texcoco.binary <- subset_samples(binary_table, Project %in% "Texcoco")
subset.texcoco.binary
sample_data(subset.texcoco.binary)

#Remove OTUs that are not present in any of the samples of the subset (OTUs that are only present either in Izta or Texcoco) if we do, do it in alpha diversity and beta diversity analysis? 
taxa_sums(subset.texcoco.binary) [1:10]

subset.texcoco.binary<- prune_taxa(taxa_sums(subset.texcoco.binary) > 0, subset.texcoco.binary)
any(taxa_sums(subset.texcoco.binary) == 0)

taxa_sums(subset.texcoco.binary) [1:10]
subset.texcoco.binary



# Relative abundance of reads per species and treatments

# Phylum
p = plot_bar(subset.texcoco.alfa , "Site", fill = "Phylum") + facet_wrap(sample_Species ~ Type) 

order_site <- list("native", "mixed", "perturbated") # re-orer the sites on x axis
p$data$Site <- factor(p$data$Site, levels = order_site)
print(p)

p + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + # remove the dividing lines that separate OTUs inside the bars
  labs(y="Abundance of sequences") # change y-axis title


# Same as above but using ggplot instead of bar_plot:

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

# re-order how site and species appear in the graph
mdata_phylum$Site <- factor(mdata_phylum$Site,levels = c("native", "mixed", "perturbated"))
mdata_phylum$sample_Species <- factor(mdata_phylum$sample_Species,levels = c("Quercus", "Juniperus"))

# Now plot by Relative Abundance of Phylum by Site, Type of sample and Plant species 
ggplot(mdata_phylum, aes(x = Site, y = Abundance, fill = Phylum)) + 
  #facet_grid(time~.) +
  geom_bar(stat = "identity")  +
  # Remove x axis title, and rotate sample lables
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
    
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("Relative abundance of sequences (Phyla > 1%)\n") +
  ggtitle("Phylum relative abundance") + facet_grid(Type ~ sample_Species)


# Same as above but considering only Asco, Basidio and Glomero
subset.phylum <- subset_taxa(subset.texcoco.alfa, Phylum %in% c("p__Ascomycota", "p__Basidiomycota", "p__Glomeromycota"))

mdata_phylum <- subset.phylum %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

# checking the dataframe that we now created
head(mdata_phylum)

# re-order how site and species appear in the graph
mdata_phylum$Site <- factor(mdata_phylum$Site,levels = c("native", "mixed", "perturbated"))
mdata_phylum$sample_Species <- factor(mdata_phylum$sample_Species,levels = c("Quercus", "Juniperus"))

# Now plot by Relative Abundance of Phylum by Site, Type of sample and Plant species 
ggplot(mdata_phylum, aes(x = Site, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity")  +   #facet_grid(time~.)
  theme(axis.title.x = element_blank(),    # Remove x axis title, and rotate sample lables
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("OTU relative abundance per Phyla (> 1%)\n") +
  facet_grid(Type ~ sample_Species)




# Relative abundance of OTUs per species and treatments (same as above but using binary table)

mdata_phylum <- subset.texcoco.binary %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

# checking the dataframe that we now created
head(mdata_phylum)

# re-order how site and species appear in the graph
mdata_phylum$Site <- factor(mdata_phylum$Site,levels = c("native", "mixed", "perturbated"))
mdata_phylum$sample_Species <- factor(mdata_phylum$sample_Species,levels = c("Quercus", "Juniperus"))

# Now plot by Relative Abundance of Phylum by Site, Type of sample and Plant species 
ggplot(mdata_phylum, aes(x = Site, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity")  +   #facet_grid(time~.)
  theme(axis.title.x = element_blank(),    # Remove x axis title, and rotate sample lables
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("OTU relative abundance per Phyla (> 1%)\n") +
  facet_grid(Type ~ sample_Species)


# Same as above but considering only Ascomycota, Basidiomycota and Glomeromycota
subset.phylum <- subset_taxa(subset.texcoco.alfa, Phylum %in% c("p__Ascomycota", "p__Basidiomycota", "p__Glomeromycota"))

mdata_phylum <- subset.phylum %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

# checking the dataframe that we now created
head(mdata_phylum)

# re-order how site and species appear in the graph
mdata_phylum$Site <- factor(mdata_phylum$Site,levels = c("native", "mixed", "perturbated"))
mdata_phylum$sample_Species <- factor(mdata_phylum$sample_Species,levels = c("Quercus", "Juniperus"))

# Now plot by Relative Abundance of Phylum by Site, Type of sample and Plant species 
ggplot(mdata_phylum, aes(x = Site, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity")  +   #facet_grid(time~.)
  theme(axis.title.x = element_blank(),    # Remove x axis title, and rotate sample lables
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("OTU relative abundance per Phyla (> 1%)\n") +
  facet_grid(Type ~ sample_Species)


# Same but by order of fungi:
mdata_order <- subset.texcoco.binary %>%
  tax_glom(taxrank = "Order") %>%                      # agglomerate at order level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.05) %>%                         # Filter out low abundance taxa
  arrange(Order)                                       # Sort data frame alphabetically by phylum

# checking the dataframe that we now created
head(mdata_order)

# re-order how site and species appear in the graph
mdata_order$Site <- factor(mdata_order$Site,levels = c("native", "mixed", "perturbated"))
mdata_order$sample_Species <- factor(mdata_order$sample_Species,levels = c("Quercus", "Juniperus"))

# Now plot by Relative Abundance of order by Site, Type of sample and Plant species
ggplot(mdata_order, aes(x = Site, y = Abundance, fill = Order)) + 
  geom_bar(stat = "identity")  +   #facet_grid(time~.)
  theme(axis.title.x = element_blank(),    # Remove x axis title, and rotate sample lables
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("OTU relative abundance per Order (> 5%)\n") +
  facet_grid(Type ~ sample_Species)



# Same but for top 50 OTUs (binary table) 

Top50OTUs <- names(sort(taxa_sums(subset.texcoco.binary), TRUE)[1:50])
ent50 <- prune_taxa(Top50OTUs, subset.texcoco.binary)
ent50


p=plot_bar(ent50, "Site", fill = "Family") +
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") + # remove the dividing lines that separate OTUs inside the bars
  labs(y="Abundance of sequences (Top 50 OTUs)") # change y-axis title

order_site <- list("native", "mixed", "perturbated") # re-orer the sites on x axis
p$data$Site <- factor(p$data$Site, levels = order_site)
print(p)


plot_bar(ent50, "sample_Species", fill = "Family", facet_grid = "Site") +
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") + # remove the dividing lines that separate OTUs inside the bars
  labs(y="Abundance of sequences (Top 50 OTUs)") # change y-axis title



mdata_phylum <- ent50 %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Family)                                      # Sort data frame alphabetically by phylum

# checking the dataframe that we now created
head(mdata_phylum)

# re-order how site and species appear in the graph
mdata_phylum$Site <- factor(mdata_phylum$Site,levels = c("native", "mixed", "perturbated"))
mdata_phylum$sample_Species <- factor(mdata_phylum$sample_Species,levels = c("Quercus", "Juniperus"))

# Now plot by Relative Abundance of Phylum by Site, Type of sample and Plant species 
ggplot(mdata_phylum, aes(x = Site, y = Abundance, fill = Family)) + 
  #facet_grid(time~.) +
  geom_bar(stat = "identity")  +
  # Remove x axis title, and rotate sample lables
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("OTU relative abundance per Family (Top 50 OTUs)") +
  facet_grid(Type ~ sample_Species)



# Same as above (Family) but for each sample
plot_bar(ent50, fill = "Family")
plot_bar(ent50, fill = "Family") + facet_wrap(~Site, scales="free_x", nrow=1)
plot_bar(ent50, fill = "Family") + facet_wrap(~sample_Species+ Site + Type, scales="free_x", nrow=1)




######################################## Subset: Izta 

# Subset data for Izta using relative abundance 
subset.izta.alfa <- subset_samples(phyloseq.rel, Site%in%c("Cueva", "Joya", "Base"))
subset.izta.alfa
sample_data(subset.izta.alfa)

sample_sums(subset.izta.alfa) [1:10] 

# Should we remove OTUs that are not present in any of the samples of the subset? (OTUs that are only present either in Izta or Texcoco) 
taxa_sums(subset.izta.alfa) [1:10]

# In case of removing OTUs not present in samples of subset:

subset.izta.alfa <- prune_taxa(taxa_sums(subset.izta.alfa) > 0, subset.izta.alfa)
any(taxa_sums(subset.izta.alfa) == 0)

taxa_sums(subset.izta.alfa) [1:10]
subset.izta.alfa

# get_taxa is returning all the OTU abundances from one sample, while get_sample is returning the abundances from all samples for one OTU.
get_taxa(subset.izta.alfa, sample_names(subset.izta.alfa)[5])[1:10]
get_sample(subset.izta.alfa, taxa_names(subset.izta.alfa)[5])[1:10]


# Relative abundance of reads per species and treatments
# melt to long format (for ggploting) 
# prune out phyla below 1% in each sample
# selecting the taxa at the level: Phylum

mdata_phylum <- subset.izta.alfa %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

# checking the dataframe that we now created
head(mdata_phylum)

# Now plot by Relative Abundance of Phylum by Site, Type of sample and Plant species 
ggplot(mdata_phylum, aes(x = Site, y = Abundance, fill = Phylum)) + 
  #facet_grid(time~.) +
  geom_bar(stat = "identity")  +
  # Remove x axis title, and rotate sample lables
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("Relative Abundance (Phyla > 1%)\n") +
  ggtitle("Phylum Relative Abundance") + facet_grid(Type ~ sample_Species)


# Subset data for izta using binary table 

subset.izta.binary <- subset_samples(binary_table, Site%in%c("Cueva", "Joya", "Base"))
subset.izta.binary
sample_data(subset.izta.binary)

#Should we remove OTUs that are not present in any of the samples of the subset? (OTUs that are only present either in Izta or Texcoco) if we do, do it in alpha diversity and beta diversity analysis? 
taxa_sums(subset.izta.binary) [1:10]

#In case of removing OTUs not present in samples of subset:

subset.izta.binary<- prune_taxa(taxa_sums(subset.izta.binary) > 0, subset.izta.binary)

subset.izta.binary

taxa_sums(subset.izta.binary) [1:10]


# Relative abundance of OTUs per species and treatments (same as above but using binary table subset )

mdata_phylum <- subset.izta.binary %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

# checking the dataframe that we now created
head(mdata_phylum)

# Now plot by Relative Abundance of Phylum by Site, Type of sample and Plant species 
ggplot(mdata_phylum, aes(x = Site, y = Abundance, fill = Phylum)) + 
  #facet_grid(time~.) +
  geom_bar(stat = "identity")  +
  # Remove x axis title, and rotate sample lables
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("Relative Abundance Binary Table (Phyla > 1%)\n") +
  ggtitle("Phylum Relative Abundance") + facet_grid(Type ~ sample_Species)


# Top 50 OTUs for relative abundance of izta 

Top50OTUs <- names(sort(taxa_sums(subset.izta.alfa), TRUE)[1:50])
ent50 <- prune_taxa(Top50OTUs, subset.izta.alfa)
plot_bar(ent50, "Site", fill = "Family")
plot_bar(ent50, "sample_Species", fill = "Family", facet_grid = "Site")

# only Ascomycota, Basidiomycota and Glomeromycota
subset.phylum <- subset_taxa(subset.izta.alfa, Phylum %in% c("p__Ascomycota", "p__Basidiomycota", "p__Glomeromycota"))
plot_bar(subset.phylum, "Site", fill = "Phylum", facet_grid = Project ~ sample_Species)

# only Glomeromycota 
subset.phylum.glomero <- subset_taxa(subset.izta.alfa, Phylum =="p__Glomeromycota")
plot_bar(subset.phylum.glomero, "Site", fill = "Phylum", facet_grid = Project ~ sample_Species)

# Phylum
plot_bar(subset.izta.alfa , "Site", fill = "Phylum") + facet_wrap(sample_Species ~ Type) 

# Family
plot_bar(ent50, fill = "Family")
plot_bar(ent50, fill = "Family") + facet_wrap(~Site, scales="free_x", nrow=1)
plot_bar(ent50, fill = "Family") + facet_wrap(~sample_Species+ Site + Type, scales="free_x", nrow=1)

# Top 50 OTUs for izta binary table 

Top50OTUs <- names(sort(taxa_sums(subset.izta.binary), TRUE)[1:50])
ent50 <- prune_taxa(Top50OTUs, subset.izta.binary)
plot_bar(ent50, "Site", fill = "Family")
plot_bar(ent50, "sample_Species", fill = "Family", facet_grid = "Site")

# only Ascomycota, Basidiomycota and Glomeromycota
subset.phylum <- subset_taxa(subset.izta.binary, Phylum %in% c("p__Ascomycota", "p__Basidiomycota", "p__Glomeromycota"))
plot_bar(subset.phylum, "Site", fill = "Phylum", facet_grid = Project ~ sample_Species)

# only Glomeromycota 
subset.phylum.glomero <- subset_taxa(subset.izta.binary, Phylum =="p__Glomeromycota")
plot_bar(subset.phylum.glomero, "Site", fill = "Phylum", facet_grid = Project ~ sample_Species)

# Phylum
plot_bar(subset.izta.binary , "Site", fill = "Phylum") + facet_wrap(sample_Species ~ Type) 

# Family
plot_bar(ent50, fill = "Family")
plot_bar(ent50, fill = "Family") + facet_wrap(~Site, scales="free_x", nrow=1)
plot_bar(ent50, fill = "Family") + facet_wrap(~sample_Species+ Site + Type, scales="free_x", nrow=1)






############################################ Alfa diversity 

# Texcoco 

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

# Observed
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


## Statistical tests 

# Create a table that gathers diversity measures
texcocodiversity<-estimate_richness(subset.texcoco.alfa, measures=c("Observed", "Fisher", "Shannon"))
data <- cbind(sample_data(subset.texcoco.alfa), texcocodiversity) #combine metadata with alpha diversity
data

# Differences between species
subset.texcoco.alfa.anova <- aov(Observed ~ Species, data =data)
summary(subset.texcoco.alfa.anova)
TukeyHSD(subset.texcoco.alfa.anova)
boxplot(Observed~Species, data = data)

# Differences between sites
subset.texcoco.alfa.anova <- aov(Observed ~ Site, data =data)
summary(subset.texcoco.alfa.anova)
TukeyHSD(subset.texcoco.alfa.anova)
boxplot(Observed~Site, data = data)

# Differences between sample types
subset.texcoco.alfa.anova <- aov(Observed ~ Type, data =data)
summary(subset.texcoco.alfa.anova)
TukeyHSD(subset.texcoco.alfa.anova)
boxplot(Observed~Type, data = data)


#two-ways

#Site/Species for Observed
subset.texcoco.alfa.anova<-aov(Observed~Species*Site, data = data)
summary(subset.texcoco.alfa.anova)
TukeyHSD(subset.texcoco.alfa.anova)
boxplot(Observed~Species*Site, data = data)

#Site/Species for Shannon (this one shows significant differences between Plant Species)
subset.texcoco.alfa.anova<-aov(Shannon~Species*Site, data = data)
summary(subset.texcoco.alfa.anova)
TukeyHSD(subset.texcoco.alfa.anova)
boxplot(Shannon~Species+Site, data = data)

#Site/Type of sample for Observed (this one shows significant differences between Type of sample)
subset.texcoco.alfa.anova <- aov(Observed ~ Site*Type, data =data)
summary(subset.texcoco.alfa.anova)
TukeyHSD(subset.texcoco.alfa.anova)
boxplot(Observed~Type*Site, data = data)

#Site/Type of sample for Shanon (this one shows significant differences between Type of sample)
subset.texcoco.alfa.anova <- aov(Shannon ~ Site*Type, data =data)
summary(subset.texcoco.alfa.anova)
TukeyHSD(subset.texcoco.alfa.anova)
boxplot(Shannon~Type*Site, data = data)


# Subset by species
subset.species <- subset_samples(subset.texcoco.alfa, Species %in% "Quercus")
subset.species

subset.species <- prune_taxa(taxa_sums(subset.species) > 0, subset.species) ## remove OTU with 0 reads
any(taxa_sums(subset.species) == 0)
subset.species

texcocodiversity<-estimate_richness(subset.species, measures=c("Observed", "Fisher", "Shannon"))
data <- cbind(sample_data(subset.species), texcocodiversity) #combine metadata with alpha diversity
data

#Site for observed 
subset.texcoco.alfa.anova <- aov(Observed ~ Site, data =data)
summary(subset.texcoco.alfa.anova)
TukeyHSD(subset.texcoco.alfa.anova)
boxplot(Observed~Site, data = data)

#Type for observed 
subset.texcoco.alfa.anova <- aov(Observed ~ Type, data =data)
summary(subset.texcoco.alfa.anova)
TukeyHSD(subset.texcoco.alfa.anova)
boxplot(Observed~Type, data = data)

#Site/Type for observed 
subset.texcoco.alfa.anova<-aov(Observed~Site*Type, data = data)
summary(subset.texcoco.alfa.anova)
TukeyHSD(subset.texcoco.alfa.anova)
boxplot(Observed~Site+Type, data = data)




############################ Izta 

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
plot_richness(subset.izta.alfa,x="Site", color = "Species", measures=("Observed"))  + geom_boxplot() +facet_wrap(~Type)

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



