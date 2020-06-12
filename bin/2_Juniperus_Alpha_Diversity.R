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



############# Exploration for direct quantitative observation/comparison of abundances 

############ Subset by project: Texcoco 

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

sum(taxa_sums(subset.texcoco.alfa))

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



# Read and OTU counts by site and by host (check also specifically for am, ecm, nomyc)


sum(taxa_sums(subset.texcoco.alfa)) # total number of reads: 3602634

sum(taxa_sums(subset.texcoco.binary)) # total number of reads: 5518

#QUERCUS ABUNDANCE ALL 

subset.quercus

sum(taxa_sums(subset.quercus)) # total number of reads: 1801292

subsetquercusnat<- subset_samples(subset.quercus, Site %in% c("native"))

sum(taxa_sums(subsetquercusnat)) # total number of reads: 900535

subsetquercusmix<- subset_samples(subset.quercus, Site %in% c("mixed"))

sum(taxa_sums(subsetquercusmix)) # total number of reads: 900757


#QUERCUS PRESENCE/ABSENCE ALL 

subset.quercus<- subset_samples(subset.texcoco.binary, Host %in% c("Quercus"))

sum(taxa_sums(subset.quercus)) # 2759

subsetquercusnat<- subset_samples(subset.quercus, Site %in% c("native"))

sum(taxa_sums(subsetquercusnat)) # 1098

subsetquercusmix<- subset_samples(subset.quercus, Site %in% c("mixed"))

sum(taxa_sums(subsetquercusmix)) # : 1661


#JUNIPERUS ABUNDANCE ALL

subset.juniperus <- subset_samples(subset.texcoco.alfa, Host %in% c("Juniperus"))

sum(taxa_sums(subset.juniperus)) # total number of reads: 1801342

subsetjuniperuspertu<- subset_samples(subset.juniperus, Site %in% c("perturbated"))

sum(taxa_sums(subsetjuniperuspertu)) # total number of reads: 900614

subsetjuniperusmix<- subset_samples(subset.juniperus, Site %in% c("mixed"))

sum(taxa_sums(subsetjuniperusmix)) # total number of reads: 900728


#JUNIPERUS PRESENCE/ABSENCE ALL

subset.juniperus <- subset_samples(subset.texcoco.binary, Host %in% c("Juniperus"))

sum(taxa_sums(subset.juniperus)) #  2759

subsetjuniperuspertu<- subset_samples(subset.juniperus, Site %in% c("perturbated"))

sum(taxa_sums(subsetjuniperuspertu)) # 1284

subsetjuniperusmix<- subset_samples(subset.juniperus, Site %in% c("mixed"))

sum(taxa_sums(subsetjuniperusmix)) # 1475


# ABUNDANCE ECM

subset.ECM <- subset_taxa(subset.texcoco.alfa, Trophic %in% c("a__ecm"))

sum(taxa_sums(subset.ECM)) # total number of reads: 1378924

# PRESENCE/ABSENCE ECM

subset.ECM <- subset_taxa(subset.texcoco.binary, Trophic %in% c("a__ecm"))

sum(taxa_sums(subset.ECM)) # 493

# ABUNDANCE AM

subset.AM <- subset_taxa(subset.texcoco.alfa, Trophic %in% c("a__am"))

sum(taxa_sums(subset.AM)) # total number of reads: 29791

# PRESENCE/ABSENCE AM

subset.AM <- subset_taxa(subset.texcoco.binary, Trophic %in% c("a__am"))

sum(taxa_sums(subset.AM)) # 272

# ABUNDANCE NOMYC

subset.NOMYC <- subset_taxa(subset.texcoco.alfa, Myc %in% c("t__nomyc"))

sum(taxa_sums(subset.NOMYC)) # total number of reads: 1374747

# PRESENCE/ABSENCE NOMYC

subset.NOMYC <- subset_taxa(subset.texcoco.binary, Myc %in% c("t__nomyc"))

sum(taxa_sums(subset.NOMYC)) # 2775

# ABUNDANCE SAP

subset.SAP <- subset_taxa(subset.texcoco.alfa, Trophic %in% c("a__sap"))

sum(taxa_sums(subset.SAP)) # total number of reads: 495061

# PRESENCE/ABSENCE SAP

subset.SAP <- subset_taxa(subset.texcoco.binary, Trophic %in% c("a__sap"))

sum(taxa_sums(subset.SAP)) # 1397



#QUERCUS ABUNDANCE ECM
#QUERCUS PRESENCE/ABSENCE ECM
#QUERCUS ABUNDANCE AM
#QUERCUS PRESENCE/ABSENCE AM
#QUERCUS ABUNDANCE NOMYC
#QUERCUS PRESENCE/ABSENCE NOMYC


#JUNIPERUS ABUNDANCE ECM
#JUNIPERUS PRESENCE/ABSENCE ECM
#JUNIPERUS ABUNDANCE AM
#JUNIPERUS PRESENCE/ABSENCE AM
#JUNIPERUS ABUNDANCE NOMYC
#JUNIPERUS PRESENCE/ABSENCE NOMYC



sample_sums(phyloseq_tables)
sample_sums(phyloseq_tables_cleaned)

tax_glom()

#ACCUMULATION CURVES / Rarefaction






# Relative abundance of reads by plant host and treatments

# Phylum
p = plot_bar(subset.texcoco.alfa , "Site", fill = "Phylum") + facet_wrap(Host ~ Type)

order_site <- list("native", "mixed", "perturbated") # re-orer the sites on x axis
p$data$Site <- factor(p$data$Site, levels = order_site)
print(p)

p + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + # remove the dividing lines that separate OTUs inside the bars
  labs(y="Relative abundance of sequences") # change y-axis title


# Same as above but using ggplot instead of bar_plot: # Relative abundance of reads by plant host and treatments

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
mdata_phylum$Host <- factor(mdata_phylum$Host,levels = c("Quercus", "Juniperus"))

# Now plot: 

ggplot(mdata_phylum, aes(x = Site, y = Abundance, fill = Phylum)) + 
  #facet_grid(time~.) +
  geom_bar(stat = "identity")  +
  # Remove x axis title, and rotate sample lables
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
    
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("Relative abundance of sequences (Phyla > 1%)\n") +
  ggtitle("Phylum relative abundance") + facet_grid(Type ~ Host)

# Plot obtained: Relative Abundance (using sequence reads) of fungi Phylum by Site, Type of sample and Host plant 


# Relative abundance of OTUs (binary table) by plant host per treatment (same as above but using binary table)

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
mdata_phylum$Host <- factor(mdata_phylum$Host,levels = c("Quercus", "Juniperus"))

# Now plot by Relative Abundance of Phylum by Site, Type of sample and Plant species 
ggplot(mdata_phylum, aes(x = Site, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity")  +   #facet_grid(time~.)
  theme(axis.title.x = element_blank(),    # Remove x axis title, and rotate sample lables
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("OTU relative abundance (> 1%)\n") +
  facet_grid(Type ~ Host)

# Plot obtained: OTU relative abundance of fungi Phylum by Site, Type of sample and Host plant 


# Repeat but using order: (Relative abundance of OTUs (binary table) by plant host and treatment)

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
mdata_order$Host <- factor(mdata_order$Host,levels = c("Quercus", "Juniperus"))

# Now plot: 
ggplot(mdata_order, aes(x = Site, y = Abundance, fill = Order)) + 
  geom_bar(stat = "identity")  +   #facet_grid(time~.)
  theme(axis.title.x = element_blank(),    # Remove x axis title, and rotate sample lables
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("OTU relative abundance (> 5%)\n") +
  facet_grid(Type ~ Host)

# Plot obtained: OTU Relative Abundance of fungal Orders by Site, Type of sample and Host plant 

# !This previous plot has something. It is not showing relative abundance? some sites show more otu abundances than others

# Abundance of TOP 50 OTUs (using binary table)

Top50OTUs <- names(sort(taxa_sums(subset.texcoco.binary), TRUE)[1:50])
ent50 <- prune_taxa(Top50OTUs, subset.texcoco.binary)
ent50

# by site
p=plot_bar(ent50, "Site", fill = "Family") +
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") + # remove the dividing lines that separate OTUs inside the bars
  labs(y="Abundance of Top 50 OTUs") # change y-axis title

order_site <- list("native", "mixed", "perturbated") # re-orer the sites on x axis
p$data$Site <- factor(p$data$Site, levels = order_site)
print(p) # Plot obtained: Abundance of TOP 50 OTUs by Family per Site 

# by host plant
p=plot_bar(ent50, "Host", fill = "Family", facet_grid = "Site") +
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") + # remove the dividing lines that separate OTUs inside the bars
  labs(y="Abundance of Top 50 OTUs") # change y-axis title

print(p) # Plot obtained: Abundance of TOP 50 OTUs by family per host plant and site 


# Relative abundance of TOP 50 OTUs by Family per host and site:

mdata_phylum <- ent50 %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Family)                                      # Sort data frame alphabetically by phylum

# checking the dataframe that we now created
head(mdata_phylum)

# re-order how site and species appear in the graph
mdata_phylum$Site <- factor(mdata_phylum$Site,levels = c("native", "mixed", "perturbated"))
mdata_phylum$Host <- factor(mdata_phylum$Host,levels = c("Quercus", "Juniperus"))

# Now plot:  
ggplot(mdata_phylum, aes(x = Site, y = Abundance, fill = Family)) + 
  #facet_grid(time~.) +
  geom_bar(stat = "identity")  +
  # Remove x axis title, and rotate sample lables
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("Top 50 OTUs relative abundance") +
  facet_grid(Type ~ Host)

# Plot obtained: Relative Abundance of TOP 50 OTUs (binary table) using Phylum by site, type and host plant 


# Abundance of OTUs (binary) by Family in TOP 50 in each sample:

plot_bar(ent50, fill = "Family") # abundance of TOP 50 OTUs at family level by site showing each sample 
plot_bar(ent50, fill = "Family") + facet_wrap(~Site, scales="free_x", nrow=1) # abundance of TOP 50 OTUs at family level by site showing each sample 




### Exploration of Texcoco data using trophic modes 

# !!! to do: remove prefix before trophic mode for final figures (with all taxonomic ranks as well)

# Abundance of reads per plant host and treatments

p = plot_bar(subset.texcoco.alfa , "Site", fill = "Trophic") + facet_wrap(Host ~ Type) 

order_site <- list("native", "mixed", "perturbated") # re-orer the sites on x axis
p$data$Site <- factor(p$data$Site, levels = order_site)
print(p)

p + geom_bar(aes(color=Trophic, fill=Trophic), stat="identity", position="stack") + # remove the dividing lines that separate OTUs inside the bars
  labs(y="Relative abundance of sequences") # change y-axis title

# Plot obtained: Abundance of sequences using trophic modes by host plant, type of sample and site


# Relative abundance of reads by plant host and treatment 

# melt to long format (for ggploting) 
# prune out phyla below 1% in each sample
# selecting the taxa at the level: Phylum

mdata_trophic <- subset.texcoco.alfa %>%
  tax_glom(taxrank = "Trophic") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Trophic)                                      # Sort data frame alphabetically by phylum

# checking the dataframe that we now created
head(mdata_trophic)

# re-order how site and species appear in the graph
mdata_trophic$Site <- factor(mdata_trophic$Site,levels = c("native", "mixed", "perturbated"))
mdata_trophic$Host <- factor(mdata_trophic$Host,levels = c("Quercus", "Juniperus"))

# Now plot: 
ggplot(mdata_trophic, aes(x = Site, y = Abundance, fill = Trophic)) + 
  #facet_grid(time~.) +
  geom_bar(stat = "identity")  +
  # Remove x axis title, and rotate sample lables
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("Relative abundance of sequences (> 1%)\n") +
  ggtitle("Trophic mode relative abundance") + facet_grid(Type ~ Host)

#Plot obtained: relative abundance of sequences using trophic mode by host plant, type of sample and site 


# Repeat but considering only ECM: relative abundance of reads for ectomycorrhizal fungi by plant host, type and site

#subset for ecm 
subset.ecto <- subset_taxa(subset.texcoco.alfa, Trophic %in% c("a__ecm"))
subset.ecto

mdata_ecm <- subset.ecto %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Family)                                      # Sort data frame alphabetically by family

# checking the dataframe that we now created
head(mdata_ecm)

# re-order how site and species appear in the graph
mdata_ecm$Site <- factor(mdata_ecm$Site,levels = c("native", "mixed", "perturbated"))
mdata_ecm$Host <- factor(mdata_ecm$Host,levels = c("Quercus", "Juniperus"))

# Now plot: 
ggplot(mdata_ecm, aes(x = Site, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity")  +   #facet_grid(time~.)
  theme(axis.title.x = element_blank(),    # Remove x axis title, and rotate sample lables
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("ECM relative abundance of sequence reads (> 1%)\n") +
  ggtitle("ECM  relative abundance") + facet_grid(Type ~ Host)

#Plot obtained: relative abundance of ECM reads at Family level by site, type of sample and plant host 


# Relative abundance of OTUs per per plant host and treatments (same as above but using binary table)

mdata_trophic <- subset.texcoco.binary %>%
  tax_glom(taxrank = "Trophic") %>%                     # agglomerate at trophic level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Trophic)                                      # Sort data frame alphabetically by trophic mode

# checking the dataframe that we now created
head(mdata_trophic)

# re-order how site and species appear in the graph
mdata_trophic$Site <- factor(mdata_trophic$Site,levels = c("native", "mixed", "perturbated"))
mdata_trophic$Host <- factor(mdata_trophic$Host,levels = c("Quercus", "Juniperus"))

# Now plot by relative abundance at trophic mode level by site, type of sample and plant host 
ggplot(mdata_trophic, aes(x = Site, y = Abundance, fill = Trophic)) + 
  geom_bar(stat = "identity")  +   #facet_grid(time~.)
  theme(axis.title.x = element_blank(),    # Remove x axis title, and rotate sample lables
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("OTU relative abundance (> 1%)\n") +
  ggtitle("Trophic mode OTUs relative abundance") + facet_grid(Type ~ Host)

# Plot obtained: Relative abundance of OTUs trophic modes by plant host, site and type of sample 


# Relative abundance of myc/nomyc OTUs per per plant host and treatments (using binary table)

mdata_myc <- subset.texcoco.binary %>%
  tax_glom(taxrank = "Myc") %>%                     # agglomerate at trophic level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Myc)                                      # Sort data frame alphabetically by trophic mode

# checking the dataframe that we now created
head(mdata_myc)

# re-order how site and species appear in the graph
mdata_myc$Site <- factor(mdata_myc$Site,levels = c("native", "mixed", "perturbated"))
mdata_myc$Host <- factor(mdata_myc$Host,levels = c("Quercus", "Juniperus"))

# Now plot by relative abundance at trophic mode level by site, type of sample and plant host 
ggplot(mdata_myc, aes(x = Site, y = Abundance, fill = Myc)) + 
  geom_bar(stat = "identity")  +   #facet_grid(time~.)
  theme(axis.title.x = element_blank(),    # Remove x axis title, and rotate sample lables
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("OTU relative abundance (> 1%)\n") +
  ggtitle("OTUs relative abundance") + facet_grid(Type ~ Host)

# Plot obtained: Relative abundance of OTUs myc/nomyc cathegory by plant host, site and type of sample 


# Same as above but considering only mycorrhizal trophic modes (binary table)
subset.myc <- subset_taxa(subset.texcoco.binary, Myc %in% c("t__myc"))
subset.myc

mdata_myc <- subset.myc %>%
  tax_glom(taxrank = "Trophic") %>%                     # agglomerate at family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Trophic)                                      # Sort data frame alphabetically by family

# checking the dataframe that we now created
head(mdata_myc)

# re-order how site and species appear in the graph
mdata_myc$Site <- factor(mdata_myc$Site,levels = c("native", "mixed", "perturbated"))
mdata_myc$Host <- factor(mdata_myc$Host,levels = c("Quercus", "Juniperus"))

# Now plot: 
ggplot(mdata_myc, aes(x = Site, y = Abundance, fill = Trophic)) + 
  geom_bar(stat = "identity")  +   #facet_grid(time~.)
  theme(axis.title.x = element_blank(),    # Remove x axis title, and rotate sample lables
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("OTU relative abundance (> 1%)\n") +
  ggtitle("Mycorrhizal OTUs relative abundance") + facet_grid(Type ~ Host)

# Plot obtained: Relative abundance of mycorrhizal OTUs by plant host, site and type of sample 


# Same as above but considering only ECM (also using binary table)
subset.ecto <- subset_taxa(subset.texcoco.binary, Trophic %in% c("a__ecm"))
subset.ecto

mdata_ecm <- subset.ecto %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Family)                                      # Sort data frame alphabetically by family

# checking the dataframe that we now created
head(mdata_ecm)

# re-order how site and species appear in the graph
mdata_ecm$Site <- factor(mdata_ecm$Site,levels = c("native", "mixed", "perturbated"))
mdata_ecm$Host <- factor(mdata_ecm$Host,levels = c("Quercus", "Juniperus"))

# Now plot: 
ggplot(mdata_ecm, aes(x = Site, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity")  +   #facet_grid(time~.)
  theme(axis.title.x = element_blank(),    # Remove x axis title, and rotate sample lables
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("OTU relative abundance (> 1%)\n") +
  ggtitle("ECM OTUs relative abundance") + facet_grid(Type ~ Host)

# Plot obtained: Relative abundance of ECM OTUs at Family level by plant host, site and type of sample 
# ! this plot is failing to show relative abundance in root samples of Juniperus perturbated 


# Repeat but considering only AM (also using binary table)

subset.am <- subset_taxa(subset.texcoco.binary, Trophic %in% c("a__am"))
subset.am

mdata_am <- subset.am %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Family)                                      # Sort data frame alphabetically by family

# checking the dataframe that we now created
head(mdata_am)

# re-order how site and species appear in the graph
mdata_am$Site <- factor(mdata_am$Site,levels = c("native", "mixed", "perturbated"))
mdata_am$Host <- factor(mdata_am$Host,levels = c("Quercus", "Juniperus"))

# Now plot: 
ggplot(mdata_am, aes(x = Site, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity")  +   #facet_grid(time~.)
  theme(axis.title.x = element_blank(),    # Remove x axis title, and rotate sample lables
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("OTU relative abundance (> 1%)\n") +
  ggtitle("AM OTUs relative abundance") + facet_grid(Type ~ Host)

# Plot obtained: Relative abundance of AM OTUs at Family level by plant host, site and type of sample 
# this plot fails to show relative abundance in soil samples of Juniperus perturbated 


# Endophytes 

subset.endo <- subset_taxa(subset.texcoco.binary, Trophic %in% c("a__endo"))
subset.endo

# Saprobes

subset.sap <- subset_taxa(subset.texcoco.binary, Trophic %in% c("a__sap"))
subset.sap 

# Parasite/pathogen fungi 

subset.par <- subset_taxa(subset.texcoco.binary, Trophic %in% c("a__par"))
subset.par

# Fungi which trophic modes that remained unresolved: unknown 
subset.unknown <- subset_taxa(subset.texcoco.binary, Trophic %in% c("a__unknown"))
subset.unknown

# Other types of mycorrhizal fungi (orchid, ericoid, )
subset.other.myc <- subset_taxa(subset.texcoco.binary, Trophic %in% c("a__otro"))
subset.other.myc

# Lichens 
subset.lic <- subset_taxa(subset.texcoco.binary, Trophic %in% c("a__lic"))
subset.lic

# Trophic mode plots for TOP 50 OTUs (binary table) 

Top50OTUs <- names(sort(taxa_sums(subset.texcoco.binary), TRUE)[1:50])
ent50 <- prune_taxa(Top50OTUs, subset.texcoco.binary)
ent50

# Plot abundance of TOP 50 OTUs 
p=plot_bar(ent50, "Site", fill = "Trophic") +
  geom_bar(aes(color=Trophic, fill=Trophic), stat="identity", position="stack") + # remove the dividing lines that separate OTUs inside the bars
  labs(y="TOP 50 OTUs abundance") # change y-axis title

order_site <- list("native", "mixed", "perturbated") # re-orer the sites on x axis
p$data$Site <- factor(p$data$Site, levels = order_site)
print(p) # OTU abundance for TOP 50 for trophic modes by site 

plot_bar(ent50, "Host", fill = "Trophic", facet_grid = "Site") +
  geom_bar(aes(color=Trophic, fill=Trophic), stat="identity", position="stack") + # remove the dividing lines that separate OTUs inside the bars
  labs(y="TOP 50 OTUs abundance") # change y-axis title
# plot obtained: OTU abundance for TOP 50 for trophic modes by site and host plant 


# OTU relative abundance for TOP 50 (using binary table)

mdata_trophic <- ent50 %>%
  tax_glom(taxrank = "Trophic") %>%                     # agglomerate at trophic level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Trophic)                                      # Sort data frame alphabetically by trophic mode

# checking the dataframe that we now created
head(mdata_trophic)

# re-order how site and species appear in the graph
mdata_trophic$Site <- factor(mdata_trophic$Site,levels = c("native", "mixed", "perturbated"))
mdata_trophic$Host <- factor(mdata_trophic$Host,levels = c("Quercus", "Juniperus"))

# Now plot: OTU relative abundance trophic mode by site, type of sample and plant host
ggplot(mdata_trophic, aes(x = Site, y = Abundance, fill = Trophic)) + 
  #facet_grid(time~.) +
  geom_bar(stat = "identity")  +
  # Remove x axis title, and rotate sample lables
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("OTU relative abundance") +
  ggtitle("TOP 50 OTUs relative abundance") + facet_grid(Type ~ Host)

#Plot obtained: TOP 50 OTUs relative abundance trophic mode by site, host plant and sample type 


# Subset using Quercus
subsetquercus <- subset_samples(subset.texcoco.binary, Host %in% c("Quercus"))
p = plot_bar(subsetquercus, "Site",  fill = "Trophic") + facet_wrap ("Type")

order_site <- list("native", "mixed") # re-order the sites on x axis
p$data$Site <- factor(p$data$Site, levels = order_site)

p + geom_bar(aes(color=Trophic, fill=Trophic), stat="identity", position="stack") + # remove the dividing lines that separate OTUs inside the bars
  labs(y="OTU abundance in Quercus") # change y-axis title



# OTU relative abundance trophic mode for Quercus 

mdata_quercus <- subsetquercus %>%
  tax_glom(taxrank = "Trophic") %>%                     # agglomerate at trophic level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Trophic)                                      # Sort data frame alphabetically by trophic mode

# checking the dataframe that we now created
head(mdata_quercus)


# Now plot: OTU relative abundance trophic mode by site, type of sample and plant host
ggplot(mdata_quercus, aes(x = Site, y = Abundance, fill = Trophic)) + 
  #facet_grid(time~.) +
  geom_bar(stat = "identity")  +
  # Remove x axis title, and rotate sample lables
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("OTU relative abundance") +
  ggtitle("Quercus: OTUs relative abundance") + facet_grid(Type ~ Host)

#Plot obtained: Quercus relative abundance trophic mode by site, host plant and sample type


# OTU relative abundance mycorrhizal trophic mode for Quercus 

mdata_quercus <- subsetquercus %>%
  tax_glom(taxrank = "Myc") %>%                     # agglomerate at trophic level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Myc)                                      # Sort data frame alphabetically by trophic mode

# checking the dataframe that we now created
head(mdata_quercus)


# Now plot: OTU relative abundance trophic mode by site, type of sample and plant host
ggplot(mdata_quercus, aes(x = Site, y = Abundance, fill = Myc)) + 
  #facet_grid(time~.) +
  geom_bar(stat = "identity")  +
  # Remove x axis title, and rotate sample lables
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("OTU relative abundance") +
  ggtitle("Quercus: OTUs relative abundance") + facet_grid(Type ~ Host)




# Subset using Juniperus
subsetjun <- subset_samples(subset.texcoco.binary, Host %in% c("Juniperus"))
p = plot_bar(subsetjun, "Site",  fill = "Trophic") + facet_wrap ("Type")

order_site <- list("mixed", "perturbated") # re-order the sites on x axis
p$data$Site <- factor(p$data$Site, levels = order_site)

p + geom_bar(aes(color=Trophic, fill=Trophic), stat="identity", position="stack") + # remove the dividing lines that separate OTUs inside the bars
  labs(y="OTU abundance in Juniperus") # change y-axis title


# OTU relative abundance trophic mode for Juniperus

mdata_jun <- subsetjun %>%
  tax_glom(taxrank = "Trophic") %>%                     # agglomerate at trophic level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Trophic)                                      # Sort data frame alphabetically by trophic mode

# checking the dataframe that we now created
head(mdata_jun)


# Now plot: OTU relative abundance trophic mode by site, type of sample and plant host
ggplot(mdata_jun, aes(x = Site, y = Abundance, fill = Trophic)) + 
  #facet_grid(time~.) +
  geom_bar(stat = "identity")  +
  # Remove x axis title, and rotate sample lables
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("OTU relative abundance") +
  ggtitle("Juniperus: OTUs relative abundance") + facet_grid(Type ~ Host)

#Plot obtained: Juniperus relative abundance trophic mode by site, host plant and sample type


# OTU relative abundance mycorrhizal trophic mode for Juniperus

mdata_jun <- subsetjun %>%
  tax_glom(taxrank = "Myc") %>%                     # agglomerate at trophic level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Myc)                                      # Sort data frame alphabetically by trophic mode

# Now plot: OTU relative abundance trophic mode by site, type of sample and plant host
ggplot(mdata_jun, aes(x = Site, y = Abundance, fill = Myc)) + 
  #facet_grid(time~.) +
  geom_bar(stat = "identity")  +
  # Remove x axis title, and rotate sample lables
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("OTU relative abundance") +
  ggtitle("Juniperus: OTUs relative abundance") + facet_grid(Type ~ Host)




### ANOVAs 

# Plotting: Quercus 
# boxplots: 

subsetquercus <- subset_samples(subset.texcoco.alfa, Host %in% c("Quercus")) # 

subsetquercus <- subset_samples(subset.texcoco.binary, Host %in% c("Quercus")) # Analysis based on this one


# fungal species richness per site  
plot_richness(subsetquercus, x="Site", color = "Type", measures=c("Observed", "Fisher", "Shannon"))  + geom_boxplot()

# how to make a subset that has only mycorrhizal taxa from Quercus? i.e. myc in Myc in taxa table and Quercus in Host in sample data

# fungal species richness per site showing only "Observed": 
plot_richness(subsetquercus, x="Site", color = "Type", measures=("Observed"))  + geom_boxplot() 

# Observed
a <- plot_richness(subsetquercus, x="Site", color = "Type", measures=("Observed"))  
ab <- a + geom_boxplot(data = a$data, aes(x=Site, y=value, color=Type, alpha=0.1)) 
ab +facet_wrap(~Type)
print(ab)

# Boxplots: Juniperus

# fungal species richness per site  
plot_richness(subsetjun, x="Site", color = "Type", measures=c("Observed", "Fisher", "Shannon"))  + geom_boxplot()

# how to make a subset that has only mycorrhizal taxa from Quercus? i.e. myc in Myc in taxa table and Quercus in Host in sample data

# fungal species richness per site showing only "Observed": 
plot_richness(subsetjun, x="Site", color = "Type", measures=("Observed"))  + geom_boxplot() 



# Create a table that gathers diversity measures and subset for mycorrhizal cathegory in all hosts

subset.myc <- subset_taxa(subset.texcoco.alfa, Myc %in% c("t__myc"))
subset.myc
myc.diversity <-estimate_richness(subset.myc, measures=c("Observed", "Fisher", "Shannon"))
data <- cbind(sample_data(subset.myc), myc.diversity) #combine metadata with alpha diversity
data


# All together
subset.myc.anova <- aov(Observed ~ Host+Site+Type, data =data)
summary(subset.myc.anova)
TukeyHSD(subset.myc.anova)

#QUERCUS MYCORRHIZAL 

# Quercus table that gathers diversity measures and subset for mycorrhizal cathegory

subset.myc <- subset_taxa(subsetquercus, Myc %in% c("t__myc"))
subset.myc
myc.diversity <-estimate_richness(subset.myc, measures=c("Observed", "Fisher", "Shannon"))
data <- cbind(sample_data(subset.myc), myc.diversity) #combine metadata with alpha diversity
data

# Quercus Observed Myc ANOVA 
subset.myc.anova <- aov(Observed ~ Site+Type, data =data)
summary(subset.myc.anova)
TukeyHSD(subset.myc.anova)
boxplot(Observed~Site*Type, data = data)

# Quercus Shannon Myc ANOVA 
subset.myc.anova <- aov(Shannon ~ Site+Type, data =data)
summary(subset.myc.anova)
TukeyHSD(subset.myc.anova)
boxplot(Shannon~Site*Type, data = data)


# QUERCUS ECTOMYCORRHIZAL

subset.myc <- subset_taxa(subsetquercus, Trophic %in% c("a__ecm"))
subset.myc
myc.diversity <-estimate_richness(subset.myc, measures=c("Observed", "Fisher", "Shannon"))
data <- cbind(sample_data(subset.myc), myc.diversity) #combine metadata with alpha diversity
data


# Quercus Observed ECM ANOVA 
subset.myc.anova <- aov(Observed ~ Site*Type, data =data)
summary(subset.myc.anova)
TukeyHSD(subset.myc.anova)
boxplot(Observed~Site*Type, data = data)

# Quercus Shannon ECM ANOVA 
subset.myc.anova <- aov(Shannon ~ Site*Type, data =data)
summary(subset.myc.anova)
TukeyHSD(subset.myc.anova)
boxplot(Shannon~Site*Type, data = data)


# QUERCUS AM 

subset.myc <- subset_taxa(subsetquercus, Trophic %in% c("a__am"))
subset.myc
myc.diversity <-estimate_richness(subset.myc, measures=c("Observed", "Fisher", "Shannon"))
data <- cbind(sample_data(subset.myc), myc.diversity) #combine metadata with alpha diversity
data


# Quercus Observed AM ANOVA 
subset.myc.anova <- aov(Observed ~ Site*Type, data =data)
summary(subset.myc.anova) # this one shows significant diff with not onluy type but with site:type
TukeyHSD(subset.myc.anova)
boxplot(Observed~Site*Type, data = data)

# Quercus Shannon AM ANOVA 
subset.myc.anova <- aov(Shannon ~ Site*Type, data =data)
summary(subset.myc.anova)
TukeyHSD(subset.myc.anova)
boxplot(Shannon~Site*Type, data = data)




#JUNIPERUS MYCORRHIZAL 

subsetjun <- subset_samples(subset.texcoco.alfa, Host %in% c("Juniperus"))

subsetjun <- subset_samples(subset.texcoco.binary, Host %in% c("Juniperus")) # Analysis based on this one


# Juniperus table that gathers diversity measures and subset for mycorrhizal cathegory

subset.myc <- subset_taxa(subsetjun, Myc %in% c("t__myc"))
subset.myc
myc.diversity <-estimate_richness(subset.myc, measures=c("Observed", "Fisher", "Shannon"))
data <- cbind(sample_data(subset.myc), myc.diversity) #combine metadata with alpha diversity
data

# Juniperus Observed Myc ANOVA 
subset.myc.anova <- aov(Observed ~ Site+Type, data =data)
summary(subset.myc.anova)
TukeyHSD(subset.myc.anova)
boxplot(Observed~Site*Type, data = data)

# Juniperus Shannon Myc ANOVA 
subset.myc.anova <- aov(Shannon ~ Site+Type, data =data)
summary(subset.myc.anova)
TukeyHSD(subset.myc.anova)
boxplot(Shannon~Site*Type, data = data)


# JUNIPERUS ECTOMYCORRHIZAL

subset.myc <- subset_taxa(subsetjun, Trophic %in% c("a__ecm"))
subset.myc
myc.diversity <-estimate_richness(subset.myc, measures=c("Observed", "Fisher", "Shannon"))
data <- cbind(sample_data(subset.myc), myc.diversity) #combine metadata with alpha diversity
data


# Juniperuss Observed ECM ANOVA 
subset.myc.anova <- aov(Observed ~ Site*Type, data =data)
summary(subset.myc.anova)
TukeyHSD(subset.myc.anova)
boxplot(Observed~Site*Type, data = data)

# Juniperus Shannon ECM ANOVA 
subset.myc.anova <- aov(Shannon ~ Site*Type, data =data)
summary(subset.myc.anova)
TukeyHSD(subset.myc.anova)
boxplot(Shannon~Site*Type, data = data)


# Juniperus AM 

subset.myc <- subset_taxa(subsetjun, Trophic %in% c("a__am"))
subset.myc
myc.diversity <-estimate_richness(subset.myc, measures=c("Observed", "Fisher", "Shannon"))
data <- cbind(sample_data(subset.myc), myc.diversity) #combine metadata with alpha diversity
data


# Juniperus Observed AM ANOVA 
subset.myc.anova <- aov(Observed ~ Site*Type, data =data)
summary(subset.myc.anova) # this one shows significant diff with not onluy type but with site:type
TukeyHSD(subset.myc.anova)
boxplot(Observed~Site*Type, data = data)

# Juniperus Shannon AM ANOVA 
subset.myc.anova <- aov(Shannon ~ Site*Type, data =data)
summary(subset.myc.anova)
TukeyHSD(subset.myc.anova)
boxplot(Shannon~Site*Type, data = data)

#############################







##################### Subset by project: Izta 

# Subset data for Izta using relative abundance of sequence reads
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

# Subset data for izta using relative abundance of OTUs (binary table, presence/absence)

subset.izta.binary <- subset_samples(binary_table, Site%in%c("Cueva", "Joya", "Base"))
subset.izta.binary
sample_data(subset.izta.binary)

#Should we remove OTUs that are not present in any of the samples of the subset? (OTUs that are only present either in Izta or Texcoco) if we do, do it in alpha diversity and beta diversity analysis? 
taxa_sums(subset.izta.binary) [1:10]

#In case of removing OTUs not present in samples of subset:

subset.izta.binary<- prune_taxa(taxa_sums(subset.izta.binary) > 0, subset.izta.binary)

subset.izta.binary

taxa_sums(subset.izta.binary) [1:10]

# get_taxa is returning all the OTU abundances from one sample, while get_sample is returning the abundances from all samples for one OTU.
get_taxa(subset.izta.binary, sample_names(subset.izta.alfa)[5])[1:10]
get_sample(subset.izta.binary, taxa_names(subset.izta.alfa)[5])[1:10]

# get_taxa is returning all the OTU abundances from one sample, while get_sample is returning the abundances from all samples for one OTU.
get_taxa(subset.izta.alfa, sample_names(subset.izta.alfa)[5])[1:10]
get_sample(subset.izta.alfa, taxa_names(subset.izta.alfa)[5])[1:10]


# Relative abundance of sequence reads by host plant and treatments 

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

# Now plot: 
ggplot(mdata_phylum, aes(x = Site, y = Abundance, fill = Phylum)) + 
  #facet_grid(time~.) +
  geom_bar(stat = "identity")  +
  # Remove x axis title, and rotate sample lables
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("Relative Abundance of sequence reads (Phyla > 1%)\n") +
  ggtitle("Phylum Relative Abundance") + facet_grid(Type ~ Host)

# Plot obtained: Relative Abundance (using sequence reads) of fungi Phylum by Site, Type of sample and Host plant 

# Relative abundance of OTUs per species and treatments (same as above but using binary table subset )

mdata_phylum <- subset.izta.binary %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

# checking the dataframe that we now created
head(mdata_phylum)

# Now plot:
ggplot(mdata_phylum, aes(x = Site, y = Abundance, fill = Phylum)) + 
  #facet_grid(time~.) +
  geom_bar(stat = "identity")  +
  # Remove x axis title, and rotate sample lables
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("OTUs Relative abundance (Phyla > 1%)\n") +
  ggtitle("OTUs Relative Abundance") + facet_grid(Type ~ Host)



# Abundance for Top 50 OTUs (izta binary table)

Top50OTUs <- names(sort(taxa_sums(subset.izta.binary), TRUE)[1:50])
ent50 <- prune_taxa(Top50OTUs, subset.izta.binary)

# Plot abundance of top 50 OTUS (binary table) at Family level by site 
plot_bar(ent50, "Site", fill = "Family")
p = plot_bar(ent50, "Host", fill = "Family", facet_grid = "Site")
p + geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") + # remove the dividing lines that separate OTUs inside the bars
  labs(y="Abundance of TOP 50 OTUs") # change y-axis title


# Relative abundance for TOP 50 OTUs at Family level 

# Relative abundance of TOP 50 OTUs by Family per host and site:

mdata_phylum <- ent50 %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Family)                                      # Sort data frame alphabetically by phylum

# checking the dataframe that we now created
head(mdata_phylum)

# Now plot:  
ggplot(mdata_phylum, aes(x = Site, y = Abundance, fill = Family)) + 
  #facet_grid(time~.) +
  geom_bar(stat = "identity")  +
  # Remove x axis title, and rotate sample lables
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("Top 50 OTUs relative abundance") +
  facet_grid(Type ~ Host)

# Plot obtained: Relative Abundance of TOP 50 OTUs (binary table) using Phylum by site, type and host plant 


# Using TROPHIC MODE


#   First no subset:

p = plot_bar(subset.izta.binary, "Site", fill = "Trophic", facet_grid = Project ~ Type)
print(p)
p + geom_bar(aes(color=Trophic, fill=Trophic), stat="identity", position="stack") + # remove the dividing lines that separate OTUs inside the bars
  labs(y="Abundance of OTUs ") # change y-axis title

# subset for OTU Abundance with AM fungi 

subset.am <- subset_taxa(subset.izta.binary, Trophic =="a__am")
p = plot_bar(subset.am, "Site", fill = "Family", facet_grid = Project ~ Type)
p + geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") + # remove the dividing lines that separate OTUs inside the bars
  labs(y="Abundance of AM OTUs") # change y-axis title

# Relative abundance of reads by plant host and treatment 

# melt to long format (for ggploting) 
# prune out phyla below 1% in each sample
# selecting the taxa at the level: Phylum

mdata_trophic <- subset.izta.alfa %>%
  tax_glom(taxrank = "Trophic") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Trophic)                                      # Sort data frame alphabetically by phylum

# Now plot: 
ggplot(mdata_trophic, aes(x = Site, y = Abundance, fill = Trophic)) + 
  #facet_grid(time~.) +
  geom_bar(stat = "identity")  +
  # Remove x axis title, and rotate sample lables
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("Relative abundance of sequences (> 1%)\n") +
  ggtitle("Trophic mode relative abundance") + facet_grid(Type ~ Host)

#Plot obtained: relative abundance of sequences using trophic mode by host plant, type of sample and site 


# OTU Relative abundance  by plant host and treatment (same as above but with binary table)

# melt to long format (for ggploting) 
# prune out phyla below 1% in each sample
# selecting the taxa at the level: Phylum

mdata_trophic <- subset.izta.binary %>%
  tax_glom(taxrank = "Trophic") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Trophic)                                      # Sort data frame alphabetically by phylum

# Now plot: 
ggplot(mdata_trophic, aes(x = Site, y = Abundance, fill = Trophic)) + 
  #facet_grid(time~.) +
  geom_bar(stat = "identity")  +
  # Remove x axis title, and rotate sample lables
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("OTUs relative abundance (> 1%)\n") +
  ggtitle("Trophic mode OTUs relative abundance") + facet_grid(Type ~ Host)

#Plot obtained: OTUs relative abundance of using trophic mode by host plant, type of sample and site 


# Same as above but considering only ECM (also using binary table)
subset.ecto <- subset_taxa(subset.izta.binary, Trophic %in% c("a__ecm"))
subset.ecto

mdata_ecm <- subset.ecto %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Family)                                      # Sort data frame alphabetically by family

# Now plot: 
ggplot(mdata_ecm, aes(x = Site, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity")  +   #facet_grid(time~.)
  theme(axis.title.x = element_blank(),    # Remove x axis title, and rotate sample lables
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("OTU relative abundance (> 1%)\n") +
  ggtitle("ECM OTUs relative abundance") + facet_grid(Type ~ Host)

# Plot obtained: Relative abundance of ECM OTUs at Family level by plant host, site and type of sample 


# Same as above but considering only AM (also using binary table)
subset.am <- subset_taxa(subset.izta.binary, Trophic %in% c("a__am"))
subset.am

mdata_am <- subset.am %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Family)                                      # Sort data frame alphabetically by family

# Now plot: 
ggplot(mdata_am, aes(x = Site, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity")  +   #facet_grid(time~.)
  theme(axis.title.x = element_blank(),    # Remove x axis title, and rotate sample lables
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("OTU relative abundance (> 1%)\n") +
  ggtitle("AM OTUs relative abundance") + facet_grid(Type ~ Host)

# Plot obtained: Relative abundance of AM OTUs at Family level by plant host, site and type of sample 

# Saprobes
subset.sap <- subset_taxa(subset.izta.binary, Trophic %in% c("a__sap"))
subset.sap

# Unknown 
subset.u <- subset_taxa(subset.izta.binary, Trophic %in% c("a__unknown"))
subset.u

# Parasite/pathogen fungi 
subset.par <- subset_taxa(subset.izta.binary, Trophic %in% c("a__par"))
subset.par

#Endo
subset.endo <- subset_taxa(subset.izta.binary, Trophic %in% c("a__endo"))
subset.endo

#Otro myc
subset.other.myc <- subset_taxa(subset.izta.binary, Trophic %in% c("a__otro"))
subset.other.myc

#Lichens
subset.lic <- subset_taxa(subset.izta.binary, Trophic %in% c("a__lic"))
subset.lic


# TOP 50 OTUs
Top50OTUs <- names(sort(taxa_sums(subset.izta.binary), TRUE)[1:50])
ent50 <- prune_taxa(Top50OTUs, subset.izta.binary)

# Abundance at Trophic mode level for the top 50 OTUS for all samples 

p = plot_bar(ent50, fill = "Trophic")
p + geom_bar(aes(color=Trophic, fill=Trophic), stat="identity", position="stack") + # remove the dividing lines that separate OTUs inside the bars
  labs(y="OTU Abundance") # change y-axis title


# OTU relative abundance for TOP 50 (using binary table)

mdata_trophic <- ent50 %>%
  tax_glom(taxrank = "Trophic") %>%                     # agglomerate at trophic level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Trophic)                                      # Sort data frame alphabetically by trophic mode

# checking the dataframe that we now created
head(mdata_trophic)

# Now plot: OTU relative abundance trophic mode by site and type of sample 
ggplot(mdata_trophic, aes(x = Site, y = Abundance, fill = Trophic)) + 
  #facet_grid(time~.) +
  geom_bar(stat = "identity")  +
  # Remove x axis title, and rotate sample lables
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  
  # add labels
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +  # modifying the legend
  ylab("OTU relative abundance") +
  ggtitle("TOP 50 OTUs relative abundance") + facet_grid(Type ~ Host)

#Plot obtained: TOP 50 OTUs relative abundance trophic mode by site, host plant and sample type 

Top50OTUs <- names(sort(taxa_sums(subset.izta.binary), TRUE)[1:50])
ent50 <- prune_taxa(Top50OTUs, subset.izta.binary)



############################################ Alfa diversity 


# Project: Texcoco 

# Plotting: 
# fungal species richness per plant species 
plot_richness(subset.texcoco.alfa,x="Host", color = "Site", shape = "Type", measures=c("Observed", "Fisher", "Shannon"))  + geom_point(size=3)

# fungal species richness per site
plot_richness(subset.texcoco.alfa,x="Site", color = "Host", shape = "Type", measures=c("Observed", "Fisher", "Shannon"))  + geom_point(size=3)

# boxplots: 

# fungal species richness per plant species 
plot_richness(subset.texcoco.alfa,x="Host", color = "Site", shape = "Type", measures=c("Observed", "Fisher", "Shannon"))  + geom_boxplot()

# fungal species richness per site
plot_richness(subset.texcoco.alfa,x="Site", color = "Host", shape = "Type", measures=c("Observed", "Fisher", "Shannon"))  + geom_boxplot()

# fungal species richness per site showing only "Observed": 
plot_richness(subset.texcoco.alfa,x="Site", color = "Host", measures=("Observed"))  + geom_boxplot() +facet_wrap(~Type)

# Observed
a <- plot_richness(subset.texcoco.alfa,x="Site", color = "Host", measures=("Observed"))  
ab <- a + geom_boxplot(data = a$data, aes(x=Site, y=value, color=Host, alpha=0.1)) 
ab +facet_wrap(~Type)

# shannon:
b <- plot_richness(subset.texcoco.alfa, x = "Site", color = "Host", measures=("Shannon")) 
bc <- b + geom_boxplot(data = b$data, aes(x=Site, y=value, color=Host, alpha=0.1))
bc + facet_wrap(~Type)

# fisher:
c <- plot_richness(subset.texcoco.alfa, x = "Site", color = "Host", measures=("Fisher")) 
cc <- c + geom_boxplot(data = c$data, aes(x=Site, y=value, color=Host, alpha=0.1))
cc + facet_wrap(~Type)


## Statistical tests 

# Create a table that gathers diversity measures
texcocodiversity<-estimate_richness(subset.texcoco.alfa, measures=c("Observed", "Fisher", "Shannon"))
data <- cbind(sample_data(subset.texcoco.alfa), texcocodiversity) #combine metadata with alpha diversity
data

# Differences between species
subset.texcoco.alfa.anova <- aov(Observed ~ Host, data =data)
summary(subset.texcoco.alfa.anova)
TukeyHSD(subset.texcoco.alfa.anova)
boxplot(Observed~Host, data = data)

# Differences between sites
subset.texcoco.alfa.anova <- aov(Observed ~ Site, data =data)
summary(subset.texcoco.alfa.anova)
TukeyHSD(subset.texcoco.alfa.anova)
boxplot(Observed~Site, data = data)

# Differences between sample types
subset.texcoco.alfa.anova <- aov(Observed ~ Type, data =data)
summary(subset.texcoco.alfa.anova) # Soil differences are significant 
TukeyHSD(subset.texcoco.alfa.anova)
boxplot(Observed~Type, data = data)

# All together
subset.texcoco.alfa.anova <- aov(Observed ~ Host+Site+Type, data =data)
summary(subset.texcoco.alfa.anova)
TukeyHSD(subset.texcoco.alfa.anova)

boxplot(Observed~Type, data = data)

#Site/Species for Shannon (this one shows significant differences between Plant Species)
subset.texcoco.alfa.anova<-aov(Shannon~Host*Site, data = data)
summary(subset.texcoco.alfa.anova)
TukeyHSD(subset.texcoco.alfa.anova)
boxplot(Shannon~Host+Site, data = data)

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



#############

# Project: Izta 

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
plot_richness(subset.izta.alfa,x="Site", color = "Host", measures=("Observed"))  + geom_boxplot() +facet_wrap(~Type)

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


############# 


#Texcoco 
  # Cathegory Myc

  # Create a table that gathers diversity measures and subset for mycorrhizal cathegory

subset.myc <- subset_taxa(subset.texcoco.alfa, Myc %in% c("t__myc"))
subset.myc
myc.diversity <-estimate_richness(subset.myc, measures=c("Observed", "Fisher", "Shannon"))
data <- cbind(sample_data(subset.myc), myc.diversity) #combine metadata with alpha diversity
data

# Differences between species
subset.myc.anova <- aov(Observed ~ Host, data =data)
summary(subset.myc.anova)
TukeyHSD(subset.myc.anova)
boxplot(Observed~Host, data = data)

# Differences between sites
subset.myc.anova <- aov(Observed ~ Site, data =data)
summary(subset.myc.anova)
TukeyHSD(subset.myc.anova)
boxplot(Observed~Site, data = data)

# Differences between sample types
subset.myc.anova <- aov(Observed ~ Type, data =data)
summary(subset.myc.anova) # Soil differences are significant 
TukeyHSD(subset.myc.anova)
boxplot(Observed~Type, data = data)

# All together
subset.myc.anova <- aov(Observed ~ Host+Site+Type, data =data)
summary(subset.myc.anova)
TukeyHSD(subset.myc.anova)

#Site/Type of sample for Observed (this one shows significant differences between Type of sample)
subset.myc.anova <- aov(Observed ~ Site*Type, data =data)
summary(subset.myc.anova)
TukeyHSD(subset.myc.anova)
boxplot(Observed~Type*Site, data = data)

#Site/Species for Shannon (this one shows signf differences between Host plant)
subset.myc.anova<-aov(Shannon~Host*Site, data = data) 
summary(subset.myc.anova)
TukeyHSD(subset.myc.anova)
boxplot(Shannon~Host+Site, data = data)

#Site/Type of sample for Shanon (this one shows significant differences between Type of sample)
subset.myc.anova <- aov(Shannon ~ Site*Type, data =data)
summary(subset.myc.anova)
TukeyHSD(subset.myc.anova)
boxplot(Shannon~Type*Site, data = data)

# Cathegory Ecto

# Create a table that gathers diversity measures and subset for ECM cathegory
subset.ecto <- subset_taxa(subset.texcoco.alfa, Trophic %in% c("a__ecm"))
subset.ecto
tex.diversity.ecto<-estimate_richness(subset.ecto, measures=c("Observed", "Fisher", "Shannon"))
data <- cbind(sample_data(subset.ecto), tex.diversity.ecto) #combine metadata with alpha diversity
data

# Differences between plant species (this one shows signf differences between host plants)
subset.ecto.anova <- aov(Observed ~ Host, data =data)
summary(subset.ecto.anova)
TukeyHSD(subset.ecto.anova)
boxplot(Observed~Host, data = data)

# Differences between sites (this one shows signf differences between sites, perturbated-native which is also referent to host plant)
subset.ecto.anova <- aov(Observed ~ Site, data =data)
summary(subset.ecto.anova)
TukeyHSD(subset.ecto.anova)
boxplot(Observed~Site, data = data)

# Differences between sample types
subset.ecto.anova <- aov(Observed ~ Type, data =data)
summary(subset.ecto.anova) # Soil differences are significant 
TukeyHSD(subset.ecto.anova)
boxplot(Observed~Type, data = data)

# All together
subset.ecto.anova <- aov(Observed ~ Host+Site+Type, data =data)
summary(subset.ecto.anova)
TukeyHSD(subset.ecto.anova)

#Site/Type of sample for Observed (this one shows significant differences between Site and Type of sample)
subset.ecto.anova <- aov(Observed ~ Site*Type, data =data)
summary(subset.ecto.anova)
TukeyHSD(subset.ecto.anova)
boxplot(Observed~Type*Site, data = data)

#Site/Species for Shannon 
subset.ecto.anova<-aov(Shannon~Host*Site, data = data)
summary(subset.ecto.anova)
TukeyHSD(subset.ecto.anova)
boxplot(Shannon~Host+Site, data = data)

#Site/Type of sample for Shanon (this one shows significant differences between Type of sample)
subset.ecto.anova <- aov(Shannon ~ Site*Type, data =data)
summary(subset.ecto.anova)
TukeyHSD(subset.ecto.anova)
boxplot(Shannon~Type*Site, data = data)


################################



# Same as above but for IZTA 

# Cathegory Myc

# Create a table that gathers diversity measures and subset for mycorrhizal cathegory

subset.myc <- subset_taxa(subset.izta.alfa, Myc %in% c("t__myc"))
subset.myc
myc.diversity <-estimate_richness(subset.myc, measures=c("Observed", "Fisher", "Shannon"))
data <- cbind(sample_data(subset.myc), myc.diversity) #combine metadata with alpha diversity
data

# Differences between sites
subset.myc.anova <- aov(Observed ~ Site, data =data)
summary(subset.myc.anova)
TukeyHSD(subset.myc.anova)
boxplot(Observed~Site, data = data)

# Differences between sample types
subset.myc.anova <- aov(Observed ~ Type, data =data)
summary(subset.myc.anova) # Soil differences are significant 
TukeyHSD(subset.myc.anova)
boxplot(Observed~Type, data = data)

# All together
subset.myc.anova <- aov(Observed ~ Site+Type, data =data)
summary(subset.myc.anova)
TukeyHSD(subset.myc.anova)

#Site/Type of sample for Observed (this one shows significant differences between Type of sample)
subset.myc.anova <- aov(Observed ~ Site*Type, data =data)
summary(subset.myc.anova)
TukeyHSD(subset.myc.anova)
boxplot(Observed~Type*Site, data = data)


#Site/Type of sample for Shanon (this one shows significant differences between Type of sample)
subset.myc.anova <- aov(Shannon ~ Site*Type, data =data)
summary(subset.myc.anova)
TukeyHSD(subset.myc.anova)
boxplot(Shannon~Type*Site, data = data)


# Cathegory Ecto

# Create a table that gathers diversity measures and subset for ECM cathegory
subset.ecto <- subset_taxa(subset.izta.alfa, Trophic %in% c("a__ecm"))
subset.ecto
izta.diversity.ecto<-estimate_richness(subset.ecto, measures=c("Observed", "Shannon"))
data <- cbind(sample_data(subset.ecto), izta.diversity.ecto) #combine metadata with alpha diversity
data

# Differences between sites (this one shows signf differences between sites, perturbated-native which is also referent to host plant)
subset.ecto.anova <- aov(Observed ~ Site, data =data)
summary(subset.ecto.anova)
TukeyHSD(subset.ecto.anova)
boxplot(Observed~Site, data = data)

# Differences between sample types
subset.ecto.anova <- aov(Observed ~ Type, data =data)
summary(subset.ecto.anova) # Soil differences are significant 
TukeyHSD(subset.ecto.anova)
boxplot(Observed~Type, data = data)

# All together
subset.ecto.anova <- aov(Observed ~ Site+Type, data =data)
summary(subset.ecto.anova)
TukeyHSD(subset.ecto.anova)

#Site/Type of sample for Observed (this one shows significant differences between Site and Type of sample)
subset.ecto.anova <- aov(Observed ~ Site*Type, data =data)
summary(subset.ecto.anova)
TukeyHSD(subset.ecto.anova)
boxplot(Observed~Type*Site, data = data)

#Site/Type of sample for Shanon (this one shows significant differences between Sites and Type of sample)
subset.ecto.anova <- aov(Shannon ~ Site*Type, data =data)
summary(subset.ecto.anova)
TukeyHSD(subset.ecto.anova)
boxplot(Shannon~Type*Site, data = data)

################################

# Export tables 

subset.texcoco.alfa
subset.texcoco.binary

subset.izta.alfa
subset.izta.binary

binary_table
binary_table_OTU2

phyloseq.rel



phyloseq.rel
taxonfiltered <- tax_table(binary_table)
otusfiltered <- otu_table(binary_table)
write.csv(otus, file='otus_filtered.csv')

