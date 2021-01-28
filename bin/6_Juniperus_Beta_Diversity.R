# Part 6. Script for Illumina MiSeq2x300 Juniperus projects: Texcoco and Izta
# In R this is step 3. Beta diversity 
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
library(mvabund)


theme_set(theme_bw())


# Set working directory to source file


source("../bin/4_Filter_otu_table.R")


# Beta diversity plots and tests 

#### Subset by project ####

# Use binary table where OTUs present in only one sample were removed

subset.texcoco.binary.beta <- subset_samples(binary_table_OTU2, Project %in% "Texcoco")
subset.texcoco.binary.beta
sample_data(subset.texcoco.binary.beta)

# Remove OTUs that are not present in any of the samples of the subset
any(taxa_sums(subset.texcoco.binary.beta) == 0)
taxa_sums(subset.texcoco.binary.beta) [1:10]
subset.texcoco.binary.beta <- prune_taxa(taxa_sums(subset.texcoco.binary.beta) > 0, subset.texcoco.binary.beta)
any(taxa_sums(subset.texcoco.binary.beta) == 0)
taxa_sums(subset.texcoco.binary.beta) [1:10]
subset.texcoco.binary.beta

#### Data vizualisation ####

# Reorder Texcoco for facet_wrap
sample_data(subset.texcoco.binary.beta)$Site = factor(sample_data(subset.texcoco.binary.beta)$Site, levels=c("native","mixed","perturbated"))

# If you wish to select data or subset by any of the sample variables, check:
sample_data(subset.texcoco.binary.beta)


# If there are samples that do not contain any fungi of the category you wish to analyze
# (e.g. ecm, am, sap) you can take out those samples for this analysis using the name of the sample or
# some variable associated with that sample like the phinchID:
trophicmode = subset_samples(subset.texcoco.binary.beta, phinchID != "")
trophicmode 

#Now you can proceed to make the subset of taxonomic group that you wish: 
selectedtrophic <- subset_taxa(trophicmode, Trophic %in% c(""))
selectedtrophic 


# Table for diversity measures and sample data
myc.diversity <-estimate_richness(selectedtrophic, measures=c("Observed"))
# verify what is there 
otu_table(selectedtrophic)
any(taxa_sums(selectedtrophic) == 0)
selectedtrophic<- prune_taxa(taxa_sums(selectedtrophic) > 0, selectedtrophic)
selectedtrophic
sample_data(selectedtrophic)

# NMDS Analysis (you can change the method for distance calculation)
nmds = distance(selectedtrophic, method = "raup")
nmds
ordination = ordinate(selectedtrophic, method = "NMDS", distance = nmds)
ordination
scores(ordination)

# New facet label names for supp variable
supp.labs <- c("native", "mixed", "disturbed")
names(supp.labs) <- c("native", "mixed", "perturbated")
sample_data(subset.texcoco.binary.beta)$Site = factor(sample_data(subset.texcoco.binary.beta)$Site, levels=c("native","mixed","disturbed"))

#Plot nmds 
p1 <- plot_ordination(selectedtrophic, ordination, color="Host", shape = "Type", title = "All fungi") + theme(aspect.ratio=1)+geom_point(size=3) 
print(p1)
p1 + facet_wrap(~Site)

#### Tests #### 

# Now that you have the plot is necessary to make some statistical tests about the community 

# PERMANOVA: Adonis
sampledf <- data.frame(sample_data(selectedtrophic)) # transform the sample information into a dataframe 

# now you can do the test with any of the variables you want to try 
adonis2( nmds ~ Host+Site, data = sampledf) 
adonis2( nmds ~ Type, data = sampledf) 
adonis2( nmds ~ Site, data = sampledf)
adonis2( nmds ~ Host, data = sampledf)
adonis2( nmds ~ Site+Type, data = sampledf)

# Now for another test you can do multivariate and univariate analysis of the composition 
# in order to do that you can group your data in different ways: 

# One way you can do the tests is by agglomerating OTUs at Family level

subset.fam<- subset_taxa(subset.texcoco.binary.beta, Trophic %in% "a__ecm") #group
subset.fam # verify
tax_table(subset.fam)

subset.fam2 <- tax_glom(subset.fam, taxrank = "Family") #agglomerate at taxonomic level
subset.fam2 #verify 
tax_table(subset.fam2)

taxa_sums(subset.fam2) [1:10] #verify what is there after subsetting, notice that otus will have a greater number when agglomerating for a higher taxonomic category 
any(taxa_sums(subset.fam2) == 0)
subset.fam2 <- prune_taxa(taxa_sums(subset.fam2) > 0, subset.fam2)
any(taxa_sums(subset.fam2) == 0)
subset.fam2

#In case you need to get top 50 (for ecm it comes out to only 31 families, for am only 14, so is not needed)
Top <- names(sort(taxa_sums(subset.fam2), TRUE)[1:50])
ent <- prune_taxa(Top, subset.fam2)
taxa_sums(ent)
taxa_names(ent)
tax_table(ent)
ent

# You can also have a subset that you can use for testing most species-rich families

subset.fam
tax_table(subset.fam)

spprich<-subset.fam %>% 
  psmelt %>%
  group_by(Family) %>% 
  summarize(OTUn = (unique(OTU) %>% length)) %>% 
  arrange(desc(OTUn))

# Now that you know which are the most species-rich families you can look at one or a subset of those you are interested in
subset.fam3<- subset_taxa(subset.fam, Family %in% c("f__ Thelephoraceae"))
subset.fam3
tax_table(subset.fam3)

# You can also test the top 50 OTUs

subset <- subset.texcoco.binary.beta
subset

taxa_sums(subset) [1:10]
any(taxa_sums(subset) == 0)
ent <- prune_taxa(taxa_sums(subset) > 0, subset)
any(taxa_sums(ent) == 0)
ent

Top <- names(sort(taxa_sums(ent), TRUE)[1:50])
ent50 <- prune_taxa(Top, ent)
ent50

taxa_names(ent50)
tax_table(ent50)
ntaxa(ent50)
taxa_sums(ent50)


# Once you have the data subsetted in whatever way you wish, you can proceed with the multivariate or univariate analysis:

# MVABUND : multivariate and univariate tests 

# check data
ent50
tax_table(ent50)
sample_data(ent50)

abun_table <- otu_table(ent50) # you need an abundance table 
abun_table

#transform to binary table in case using agglomeration (e.g. at family level using tax_glom)
abun_table = transform_sample_counts(abun_table, function(x, minthreshold=0){
  x[x > minthreshold] <- 1
  return(x)})
head(otu_table(abun_table))

#if not doing previous line, procede:

abun_table = t (abun_table) 
meta_table <- sample_data(ent50)  # you need the samples information
meta_table

abun_sample_table <- data.frame(abun_table, meta_table) #convert it all into a data frame
abun_sample_table

mvabund_table <- mvabund(abun_sample_table [,1:50]) #format table for mvabund. The number must be equal to the number of OTUs/families/... you are using
mvabund_table

mod2 <- manyglm(mvabund_table ~ abun_sample_table$Site , family="negative_binomial") #now fit to the model
plot(mod2)

#manyglm <- manyglm(mvabund ~ pH * factor(Elevation), family="binomial")

anova(mod2) # now the stats
anova(mod2, p.uni="adjusted") # now the stats for the univariate 




