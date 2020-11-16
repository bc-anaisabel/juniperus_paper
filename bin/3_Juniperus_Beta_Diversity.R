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
library(mvabund)
library(randomForest)
library(knitr)
library(gplots)
library(VennDiagram)
library(bipartite)
library(sna)
library(igraph)
library(ggnetwork)
library(statnet.common)
library(network)
library(pals)

theme_set(theme_bw())


# Set working directory to source file


source("../bin/1_Filter_otu_table.R")


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

# Select data or subset 
sample_data(subset.texcoco.binary.beta)


# Take out samples with no fungi present for ECM and for AM 
trophicmode = subset_samples(subset.texcoco.binary.beta, phinchID != "H8-J1r-2-2621902")
trophicmode 

#Subset if necessary
selectedtrophic <- subset_taxa(trophicmode, Trophic %in% c("a__ecm"))
selectedtrophic 


# Table for diversity measures and sample data
myc.diversity <-estimate_richness(selectedtrophic, measures=c("Observed"))

otu_table(selectedtrophic)
any(taxa_sums(selectedtrophic) == 0)
selectedtrophic<- prune_taxa(taxa_sums(selectedtrophic) > 0, selectedtrophic)
selectedtrophic
sample_data(selectedtrophic)

# NMDS Texcoco (change method for distance)
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
p1 <- plot_ordination(selectedtrophic, ordination, color="Host", shape = "Type", title = "ECM fungi") + theme(aspect.ratio=1)+geom_point(size=3) 
print(p1)
p1 + facet_wrap(~Site)

#### Tests #### 

#PERMANOVA: Adonis
sampledf <- data.frame(sample_data(selectedtrophic))

adonis2( nmds ~ Host+Site, data = sampledf) 
adonis2( nmds ~ Type, data = sampledf) 
adonis2( nmds ~ Site, data = sampledf)
adonis2( nmds ~ Host, data = sampledf)
adonis2( nmds ~ Site+Type, data = sampledf)

#MVABUND : multivariate and univariate tests 

# Agglomerating OTUs at Family level

subset.fam<- subset_taxa(subset.texcoco.binary.beta, Trophic %in% "a__sap") 
subset.fam
tax_table(subset.fam)

subset.fam <- tax_glom(subset.fam, taxrank = "Family")
subset.fam
tax_table(subset.fam)

taxa_sums(subset.fam) [1:10]
any(taxa_sums(subset.fam) == 0)
subset.fam <- prune_taxa(taxa_sums(subset.fam) > 0, subset.fam)
any(taxa_sums(subset.fam) == 0)
subset.fam

#In case you need to get top 50 (for ecm it comes out to only 31 families, for am only 14, so is not needed)
Top <- names(sort(taxa_sums(subset.fam), TRUE)[1:50])
ent <- prune_taxa(Top, subset.fam)

taxa_sums(ent)
taxa_names(ent)
tax_table(ent)
ent


# For selecting most species-rich families


subset
tax_table(subset)

spprich<-subset %>% 
  psmelt %>%
  group_by(Family) %>% 
  summarize(OTUn = (unique(OTU) %>% length)) %>% 
  arrange(desc(OTUn))

subset<- subset_taxa(subset, Family %in% c("f__ Herpotrichiellaceae"))
subset
tax_table(subset)

# For top 50 OTUs

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


#mvabund

# check data
ent50
tax_table(ent50)
sample_data(ent50)

abun_table <- otu_table(ent50) 
abun_table

#transform to binary table in case using agglomeration (e.g. at family level using tax_glom)
abun_table = transform_sample_counts(abun_table, function(x, minthreshold=0){
  x[x > minthreshold] <- 1
  return(x)})
head(otu_table(abun_table))

#if not, procede:

abun_table = t (abun_table)

meta_table <- sample_data(ent50) 
meta_table

abun_sample_table <- data.frame(abun_table, meta_table)
abun_sample_table

mvabund_table <- mvabund(abun_sample_table [,1:73]) #format table 
mvabund_table

mod2 <- manyglm(mvabund_table ~ abun_sample_table$Site , family="negative_binomial")
plot(mod2)

#manyglm <- manyglm(mvabund ~ pH * factor(Elevation), family="binomial")

anova(mod2)
anova(mod2, p.uni="adjusted")






#### Create table frequency mycorrhizal fungi in each Host####

subset.myc <- subset_taxa(subset.texcoco.binary.beta, Trophic %in% c("a__am"))
subset.myc <- subset_samples(subset.myc, Type %in% c("root"))
subset.myc

any(taxa_sums(subset.myc) == 0)
subset.myc <- prune_taxa(taxa_sums(subset.myc) > 0, subset.myc)
subset.myc

tax_table(subset.myc)
otu_table(subset.myc)


#make a dataframe with OTUs as column and samples as rows
predictors <- t(otu_table(subset.myc))
dim(predictors)

predictors <- as.table(predictors)


predictors<- as.data.frame(predictors)
dim(predictors)

Host <- as.factor(sample_data(subset.myc)$Host)
Site <- as.factor(sample_data(subset.myc)$Site)

predictors <- data.frame(predictors,Host,Site)

tax.table<-as.data.frame(tax_table(subset.myc))
tax.table
match.id <- match(predictors$Var2, rownames(tax.table))

predictors$Family <- tax.table$Family[match.id]
predictors$Trophic <- tax.table$Trophic[match.id]
predictors$Genus <- tax.table$Genus[match.id] 

predictors <- subset(predictors, predictors$Freq == 1)
  
sort(table(predictors$Var2))

sharedotusam<-write.csv(predictors, file = "sharedotusacm.csv")

###

#Random forest 

#make a dataframe of training data with OTUs as column and samples as rows
predictors <- t(otu_table(subset.myc))
dim(predictors)

#make a column for the outcome/response variable
response <- as.factor(sample_data(subset.myc)$Site)

#combine them into 1 data frame
rf.data <- data.frame(response, predictors)

set.seed(2)
fungi.classify <- randomForest(response~., data = rf.data, ntree = 400)
print(fungi.classify)

#what variables are stored in the output?
names(fungi.classify)

#make a data frame with predictor names and their importance
imp <-importance(fungi.classify)
imp <-data.frame(predictors = rownames(imp), imp)

#Order the predictor levels by importance
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

#Select the top predictors
imp.50 <- imp.sort[1:20, ]

#ggplot

ggplot(imp.50, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") + 
  coord_flip() +
  ggtitle("Most important AM OTUS (roots) for classifying samples/n into Sites")


# What are those OTUs?
otunames <- imp.50$predictors
r <- rownames(tax_table(selectedtrophic)) %in% otunames
kable(tax_table(selectedtrophic)[r, ])









#### network using sna #### 

#Select data
subset <- subset_taxa(subset.texcoco.binary.beta, Trophic %in% c("a__am"))
subset<- subset_samples(subset, Type %in% "root")
subset

any(taxa_sums(subset) == 0)
subset <- prune_taxa(taxa_sums(subset) > 0, subset)
subset

tax_table(subset)
otu_table(subset)


#Merge by category 
nw <-merge_samples(subset, group = "Host")
network_host <- as.data.frame(otu_table(nw))
is.matrix(network_host)

# vectors
taxa_names(subset)
# color_vector <- c()
color<-as.matrix(tax_table(subset))
color<-as.data.frame(color)
color<- color$Family


#Gplot
gplot(network_host, thresh = 0.2, displaylabels = TRUE, vertex.col = color)
network_host

#Merge more than one category 

sample_variables(subset)

variable1 = as.character(get_variable(subset, "Host"))
variable2 = as.character(get_variable(subset, "Site"))

sample_data(subset)$NewPastedVar <- mapply(paste0, variable1, variable2, 
                                           collapse = "_")
nw2<- merge_samples(subset, "NewPastedVar")

network_host <- as.data.frame(otu_table(nw2))
is.matrix(network_host)



#Gplot
gplot(network_host, thresh = 0.2, displaylabels = TRUE, usearrows=FALSE, 
      legend(x=1,y=-1, color_vector, pch=21, col = "#777777", 
             pt.cex=2, cex=.8, bty="n", ncol=1), vertex.col = color$Family)


gplot(network_host, thresh = 0.2, displaylabels = TRUE, label = color$Family, usearrows=FALSE, 
      legend(x=1,y=-1, color_vector, pch=21, col = "#777777", 
             pt.cex=2, cex=.8, bty="n", ncol=1), vertex.col = color$Family)



#Other format
abc<-t(network_host)

#Plot
plotweb(network_host)

#Plot
visweb(network_host)

network_host

#Other type of plot 

pal2 <-polychrome(27)

par(mfrow=c(1,2), xpd=T)

gplot(as.one.mode(network_host),
      displaylabels = TRUE,
      gmode="graph",
      label.cex=1, vertex.col = color, vertex.cex=1)

palette(polychrome(n=27))
gplot(network_host, gmode="graph", jitter=FALSE,
      displaylabels=TRUE, 
      boxed.labels=FALSE, label.pos=1, label.cex=1, vertex.cex=2,
      vertex.col= as.numeric(color$Family))
  


#### Networks using igraph ####



