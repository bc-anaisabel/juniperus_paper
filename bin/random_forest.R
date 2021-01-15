#Random forests
#For Juniperus data 2020
# Script for Illumina MiSeq2x300 Juniperus project


#Get everything ready, libraries and calling previous scripts

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


source("../bin/4_Filter_otu_table.R")


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

#First you need to create something like a frequency table

#### Create table frequency mycorrhizal fungi in each Host####

subset.myc <- subset_taxa(subset.texcoco.binary.beta, Trophic %in% c("a__ecm"))
subset.myc <- subset_samples(subset.myc, Type %in% c("root"))
subset.myc

any(taxa_sums(subset.myc) == 0)
subset.myc <- prune_taxa(taxa_sums(subset.myc) > 0, subset.myc)
subset.myc

tax_table(subset.myc)
otu_table(subset.myc)
sample_data(subset.myc)

#make a dataframe with OTUs as column and samples as rows
predictors <- t(otu_table(subset.myc))
dim(predictors)

predictors <- as.table(predictors)


predictors<- as.data.frame(predictors)
dim(predictors)

Host <- as.factor(sample_data(subset.myc)$Host)
Site <- as.factor(sample_data(subset.myc)$Site)

predictors <- data.frame(predictors,Host,Site)

head(predictors)

predictors$HS <- paste0(predictors$Host, predictors$Site)

tax.table<-as.data.frame(tax_table(subset.myc))
tax.table
match.id <- match(predictors$Var2, rownames(tax.table))

predictors$Family <- tax.table$Family[match.id]
predictors$Trophic <- tax.table$Trophic[match.id]
predictors$Genus <- tax.table$Genus[match.id] 

predictors <- subset(predictors, predictors$Freq == 1)

sort(table(predictors$Var2))

sharedotusam<-write.csv(predictors, file = "sharedotusacm.csv")


#### Random forest #### 

#make a dataframe of training data with OTUs as column and samples as rows
predictors <- t(otu_table(subset.myc))
dim(predictors)

#predictors <- as.table(predictors)
#predictors<- as.data.frame(predictors)
#dim(predictors)

#Host <- as.factor(sample_data(subset.myc)$Host)
#Site <- as.factor(sample_data(subset.myc)$Site)

#predictors <- data.frame(predictors,Host,Site)
#head(predictors)

#predictors$HS <- paste0(predictors$Host, predictors$Site)
#head(predictors)

#make a column for the outcome/response variable
response <- as.factor(sample_data(subset.myc)$Host)

#combine them into 1 data frame
rf.data <- data.frame(response, predictors)

set.seed(2)
fungi.classify <- randomForest(response ~., data = rf.data, ntree = 500)
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
  ggtitle("Most important AM OTUS (roots) for classifying samples/n into Host/Sites")


# What are those OTUs?
otunames <- imp.50$predictors
r <- rownames(tax_table(selectedtrophic)) %in% otunames
kable(tax_table(selectedtrophic)[r, ])


# Trying random forest with four levels for host/site category

#make a dataframe of training data with OTUs as column and samples as rows
predictors <- t(otu_table(subset.myc))
dim(predictors)

predictors <- as.table(predictors)
predictors<- as.data.frame(predictors)
dim(predictors)

Host <- as.factor(sample_data(subset.myc)$Host)
Site <- as.factor(sample_data(subset.myc)$Site)

predictors <- data.frame(predictors,Host,Site)
head(predictors)

predictors$HS <- paste0(predictors$Host, predictors$Site)
head(predictors)

#make a column for the outcome/response variable
response <- as.factor(sample_data(predictors)$HS)

#combine them into 1 data frame
rf.data <- data.frame(response, predictors)

set.seed(2)
fungi.classify <- randomForest(response ~., data = rf.data, ntree = 500)
print(fungi.classify)
