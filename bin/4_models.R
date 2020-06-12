# Script for Illumina MiSeq2x300 Juniperus projects: Texcoco and Izta
# Part.4 Models
# May, 2020 

#### Set-up ####

library("phyloseq")
library("vegan")
library("ggplot2")
library("plyr"); packageVersion("plyr")
library("RColorBrewer"); packageVersion("RColorBrewer")
library("plotly"); packageVersion("plotly")
library("htmltools"); packageVersion("htmltools")
library("DT"); packageVersion("DT")
library(dplyr)
library(tibble)
library("mvabund")


# Set working directory to source file


source("../bin/1_Filter_otu_table.R")

# Tables obtained from 1_Filter_otu_table for further analysis:

# phyloseq.rel (OTU table with relative abundances - corrected in count) 
# binary_table
# phyloseq.rel_OTU2 
# binary_table_OTU2 (corrected binary table, i.e. removed OTUs that appear only in 1 sample)

juniperus.rel <- phyloseq.rel  #change name of sequence abundance table for this script 

juniperus.rel
# juniperus.rel_OTU2 # TO FIX! (Right now calculate read sum and not taxa sum)

#### Separar por proyecto binomial tables:

#Texcoco 
texcoco.binary.models <- subset_samples(binary_table_OTU2, Project %in% "Texcoco")
texcoco.binary.models
sample_data(texcoco.binary.models)

# Remove OTUs that are not present in any of the samples of the subset (OTUs that are only present either in Izta or Texcoco) 
any(taxa_sums(texcoco.binary.models) == 0)
taxa_sums(texcoco.binary.models) [1:10]

texcoco.binary.models <- prune_taxa(taxa_sums(texcoco.binary.models) > 0, texcoco.binary.models)

any(taxa_sums(texcoco.binary.models) == 0)
taxa_sums(texcoco.binary.models) [1:10]
texcoco.binary.models




#Izta
izta.binary.models <- subset_samples(binary_table_OTU2, Project %in% "Izta")
izta.binary.models
sample_data(izta.binary.models)

# Remove OTUs that are not present in any of the samples of the subset (OTUs that are only present either in Izta or Texcoco) 
any(taxa_sums(izta.binary.models) == 0)
taxa_sums(izta.binary.models) [1:10]

izta.binary.models <- prune_taxa(taxa_sums(izta.binary.models) > 0, izta.binary.models)

any(taxa_sums(izta.binary.models) == 0)
taxa_sums(izta.binary.models) [1:10]
izta.binary.models


#### Models using metadata: Texcoco only soil samples ####

soil <- subset_samples(texcoco.binary.models, Type %in% "soil") # phyloseq object subset with only soil samples
soil

binary_table <- otu_table(soil) #just otu table
binary_table

meta_table <- sample_data(soil) #only sample data and only from soil samples
meta_table

env <- c("pH","Pdis","Ca","Mg","K","Na","H","SoilM") # Al removed becuase of missing data
micro <- c("Ca","Mg","K","Na","H")

#Transform columns to numeric values, each column separately 

env_table <- subset(meta_table, select = env)
env_table <-data.frame(env_table) #use as data frame to be able to transform to numeric
env_table[, 8]  <- as.numeric(env_table[, 8])
env_table[, 7]  <- as.numeric(env_table[, 7])
env_table[, 6]  <- as.numeric(env_table[, 6])
env_table[, 5]  <- as.numeric(env_table[, 5])
env_table[, 4]  <- as.numeric(env_table[, 4])
env_table[, 3]  <- as.numeric(env_table[, 3])
env_table[, 2]  <- as.numeric(env_table[, 2])
env_table[, 1]  <- as.numeric(env_table[, 1])

### env_table[, 1:8]  <- as.numeric(env_table[, 1:8])

#Transforming to log scale 

env_table_scale = scale(log(env_table), center=TRUE,scale = TRUE) # no missing data allowed
env_table_scale
env_table_scale<- data.frame(env_table_scale) #must be transformed to data frame again 

data=data.frame(env_table_scale, meta_table$Site, meta_table$Host, row.names=rownames(meta_table))
data
colnames(data) [9] <- "Site"
colnames(data) [10] <- "Host"
data

micro_table <- subset(meta_table, select = micro)
micro_table

# Distance-based redundancy analysis (dbRDA) is an ordination method similar to Redundancy Analysis (rda),
# but it allows non-Euclidean dissimilarity indices, such as Manhattan or Bray-Curtis, Raup-Crick distance. 
# Despite this non-Euclidean feature, the analysis is strictly linear and metric. 

binary_table = t (binary_table)
colnames(binary_table)
raup <- vegdist(binary_table,method="raup")
head(raup)

# null model
raup.min<-capscale (raup ~ 1, env_table_scale)
plot(raup.min)


# model including all env variables - dbrda(raup ~ ., data) is the same!
raup.all <-dbrda (raup ~ pH + Pdis + Ca + Mg + K + Na + H + SoilM, env_table_scale)
plot(raup.all)

# hypothesis model
raup.test<-dbrda(raup ~ pH + Pdis, env_table_scale)

# which is the best model?
# if both are significant, you choose the simplest model (raup.test)
anova (raup.all, raup.test)


# OrdiR2step: To select the significant explanatory variables
# The criteria for including a variable in the model is based on adjusted variation (R2adj) and Pvalues
# explained by the selected variables and their comparison with R2adj explained by the global model (with all variables);
# if the new variable is not significant or the R2adj of the model including this new variable would exceed the R2adj of the global model,
# the selection will be stopped.

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction =  "both") # Two factors are significant: Pdis (0.046 Pvalue) and then Ca (0.068 Pvalue)
RsquareAdj (raup.min)$adj.r.squared # adjusted R2 explained by all variables


raup.db # new model with selected variables based on on adjusted variation: 
plot(raup.db)
summary(raup.db)
RsquareAdj (raup.db)$adj.r.squared
anova(raup.db)
anova(raup.db, by="term") 
anova(raup.db, by="axis")


# model including host 
raup.min<-capscale (raup ~ 1, data)
raup.all <-dbrda (raup ~ pH + Pdis + Ca + Mg + K + Na + H + SoilM + Host, data) 
  

# OrdiR2step: To select the significant explanatory variables
# The criteria for including a variable in the model is based on adjusted variation (R2adj) and Pvalues
# explained by the selected variables and their comparison with R2adj explained by the global model (with all variables);
# if the new variable is not significant or the R2adj of the model including this new variable would exceed the R2adj of the global model,
# the selection will be stopped.

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction =  "both")
RsquareAdj (raup.min)$adj.r.squared # adjusted R2 explained by all variables


raup.db # new model with selected variables based on on adjusted variation:  
plot(raup.db)
summary(raup.db)
RsquareAdj (raup.db)$adj.r.squared
anova(raup.db)
anova(raup.db, by="term") 
anova(raup.db, by="axis")

# model including site 

raup.min<-capscale (raup ~ 1, data)
raup.all <-dbrda (raup ~ pH + Pdis + Ca + Mg + K + Na + H + SoilM + Site, data) 


# OrdiR2step: To select the significant explanatory variables
# The criteria for including a variable in the model is based on adjusted variation (R2adj) and Pvalues
# explained by the selected variables and their comparison with R2adj explained by the global model (with all variables);
# if the new variable is not significant or the R2adj of the model including this new variable would exceed the R2adj of the global model,
# the selection will be stopped.

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction =  "both")
RsquareAdj (raup.min)$adj.r.squared # adjusted R2 explained by all variables


raup.db # new model with selected variables based on on adjusted variation:  
plot(raup.db)
summary(raup.db)
RsquareAdj (raup.db)$adj.r.squared
anova(raup.db)
anova(raup.db, by="term") 
anova(raup.db, by="axis")



#### Models: Texcoco using Site, Host and Type as factors ####
# (Repeat integrating also root samples and instead of using environmental variables use only treatments) 

texcoco.binary.models

binary_table <- otu_table(texcoco.binary.models) #just otu table
binary_table

meta_table <- sample_data(texcoco.binary.models) #only sample data and only from soil samples
meta_table

Host <- c("Host") 
Site <- c("Site")
Type <- c("Type")
HS <- c("Site", "Host")
HST <- c("Site", "Host", "Type")

#Transform columns to numeric values, each column separately 

env_table <- subset(meta_table, select = HST)
env_table <-data.frame(env_table) #use as data frame to be able to transform to numeric
env_table[, 3]  <- as.numeric(env_table[, 3])
env_table[, 2]  <- as.numeric(env_table[, 2])
env_table[, 1]  <- as.numeric(env_table[, 1])


#Transforming to log scale 

env_table_scale = scale(log(env_table), center=TRUE,scale = TRUE) # no missing data allowed
env_table_scale
env_table_scale<- data.frame(env_table_scale) #must be transformed to data frame again 

data=data.frame(env_table_scale, meta_table$Site, meta_table$Host, meta_table$Type, row.names=rownames(meta_table))
data
colnames(data) [4] <- "Site"
colnames(data) [5] <- "Host"
colnames(data) [6] <- "Type"
data


# Distance-based redundancy analysis (dbRDA) is an ordination method similar to Redundancy Analysis (rda),
# but it allows non-Euclidean dissimilarity indices, such as Manhattan or Bray-Curtis, Raup-Crick distance. 
# Despite this non-Euclidean feature, the analysis is strictly linear and metric. 

binary_table = t (binary_table)
colnames(binary_table)
raup <- vegdist(binary_table,method="raup")
head(raup)

# null model
raup.min<-capscale (raup ~ 1, data)
plot(raup.min)

# model including host, site and type - dbrda(raup ~ ., data) is the same!
raup.all <-dbrda (raup ~ Host + Site + Type, data)
plot(raup.all)

# OrdiR2step: To select the significant explanatory variables
# The criteria for including a variable in the model is based on adjusted variation (R2adj) and Pvalues
# explained by the selected variables and their comparison with R2adj explained by the global model (with all variables);
# if the new variable is not significant or the R2adj of the model including this new variable would exceed the R2adj of the global model,
# the selection will be stopped.

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction =  "both") 
RsquareAdj (raup.min)$adj.r.squared # adjusted R2 explained by all variables


raup.db # new model with selected variables based on on adjusted variation:  
plot(raup.db)
summary(raup.db)
RsquareAdj (raup.db)$adj.r.squared
anova(raup.db)
anova(raup.db, by="term") 
anova(raup.db, by="axis")


#### Models: Texcoco using only Quercus and type and site as factors #### 

#Only Quercus samples with all fungi 
quercus <- subset_samples(texcoco.binary.models, Host %in% "Quercus") 
quercus

binary_table <- otu_table(quercus) #just otu table
binary_table

meta_table <- sample_data(quercus) 
meta_table


Host <- c("Host") 
Site <- c("Site")
Type <- c("Type")
ST <- c("Site", "Type")

#Transform columns to numeric values, each column separately 

env_table <- subset(meta_table, select = ST)
env_table <-data.frame(env_table) #use as data frame to be able to transform to numeric
env_table[, 2]  <- as.numeric(env_table[, 2])
env_table[, 1]  <- as.numeric(env_table[, 1])

### env_table[, 1:8]  <- as.numeric(env_table[, 1:8])

#Transforming to log scale 

env_table_scale = scale(log(env_table), center=TRUE,scale = TRUE) # no missing data allowed
env_table_scale
env_table_scale<- data.frame(env_table_scale) #must be transformed to data frame again 

data=data.frame(env_table_scale, meta_table$Site, meta_table$Host, row.names=rownames(meta_table))
data
colnames(data) [3] <- "Site"
colnames(data) [4] <- "Host"
data

# Distance-based redundancy analysis (dbRDA) is an ordination method similar to Redundancy Analysis (rda),
# but it allows non-Euclidean dissimilarity indices, such as Manhattan or Bray-Curtis, Raup-Crick distance. 
# Despite this non-Euclidean feature, the analysis is strictly linear and metric. 

binary_table = t (binary_table)
colnames(binary_table)
raup <- vegdist(binary_table,method="raup")
head(raup)

# null model
raup.min<-capscale (raup ~ 1, data)
plot(raup.min)


# model including all variables - dbrda(raup ~ ., data) is the same!
raup.all <-dbrda (raup ~ Site + Type, data)
plot(raup.all)

#raup.all <-dbrda (raup ~ (Site/Type), env_table_scale)
#plot(raup.all)

# OrdiR2step: To select the significant explanatory variables
# The criteria for including a variable in the model is based on adjusted variation (R2adj) and Pvalues
# explained by the selected variables and their comparison with R2adj explained by the global model (with all variables);
# if the new variable is not significant or the R2adj of the model including this new variable would exceed the R2adj of the global model,
# the selection will be stopped.

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction =  "both") 
RsquareAdj (raup.min)$adj.r.squared # adjusted R2 explained by all variables


raup.db # new model with selected variables based on on adjusted variation: 
plot(raup.db)
summary(raup.db)
RsquareAdj (raup.db)$adj.r.squared
anova(raup.db)
anova(raup.db, by="term") 
anova(raup.db, by="axis")

mod <- varpart(raup, ~Site, ~Type, data = data)
mod
showvarparts(2)
plot(mod, cutoff = -Inf, Xnames = c('Site', 'Type'))


#### Models: Texcoco using only Juniperus and type and site as factors #### 
juniperus <- subset_samples(texcoco.binary.models, Host %in% "Juniperus") 
juniperus

binary_table <- otu_table(juniperus) #just otu table
binary_table

meta_table <- sample_data(juniperus) 
meta_table

Host <- c("Host") 
Site <- c("Site")
Type <- c("Type")
ST <- c("Site", "Type")

#Transform columns to numeric values, each column separately 

env_table <- subset(meta_table, select = ST)
env_table <-data.frame(env_table) #use as data frame to be able to transform to numeric
env_table[, 2]  <- as.numeric(env_table[, 2])
env_table[, 1]  <- as.numeric(env_table[, 1])

### env_table[, 1:8]  <- as.numeric(env_table[, 1:8])

#Transforming to log scale 

env_table_scale = scale(log(env_table), center=TRUE,scale = TRUE) # no missing data allowed
env_table_scale
env_table_scale<- data.frame(env_table_scale) #must be transformed to data frame again 

data=data.frame(env_table_scale, meta_table$Site, meta_table$Host, row.names=rownames(meta_table))
data
colnames(data) [3] <- "Site"
colnames(data) [4] <- "Host"
data

# Distance-based redundancy analysis (dbRDA) is an ordination method similar to Redundancy Analysis (rda),
# but it allows non-Euclidean dissimilarity indices, such as Manhattan or Bray-Curtis, Raup-Crick distance. 
# Despite this non-Euclidean feature, the analysis is strictly linear and metric. 

binary_table = t (binary_table)
colnames(binary_table)
raup <- vegdist(binary_table,method="raup")
head(raup)

# null model
raup.min<-capscale (raup ~ 1, data)
plot(raup.min)


# model including all variables - dbrda(raup ~ ., data) is the same!
raup.all <-dbrda (raup ~ Site + Type, data)
plot(raup.all)


# OrdiR2step: To select the significant explanatory variables
# The criteria for including a variable in the model is based on adjusted variation (R2adj) and Pvalues
# explained by the selected variables and their comparison with R2adj explained by the global model (with all variables);
# if the new variable is not significant or the R2adj of the model including this new variable would exceed the R2adj of the global model,
# the selection will be stopped.

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction =  "both") 
RsquareAdj (raup.min)$adj.r.squared # adjusted R2 explained by all variables


raup.db # new model with selected variables based on on adjusted variation: 
plot(raup.db)
summary(raup.db)
RsquareAdj (raup.db)$adj.r.squared
anova(raup.db)
anova(raup.db, by="term") 
anova(raup.db, by="axis")

mod <- varpart(raup, ~Site, ~Type, data = data)
mod
showvarparts(2)
plot(mod, cutoff = -Inf, Xnames = c('Site', 'Type'))


#### Models: Texcoco using trophic modes ####


# Only mycorrhizal samples "myc" 
myc <- subset_taxa(texcoco.binary.models, Myc %in% "t__myc") 
myc

binary_table <- otu_table(myc) #just otu table
binary_table

meta_table <- sample_data(myc) #only sample data 
meta_table

Host <- c("Host") 
Site <- c("Site")
Type <- c("Type")
HST <- c("Host", "Site", "Type")


#Transform columns to numeric values, each column separately 

env_table <- subset(meta_table, select = HST)
env_table <-data.frame(env_table) #use as data frame to be able to transform to numeric
env_table[, 3]  <- as.numeric(env_table[, 3])
env_table[, 2]  <- as.numeric(env_table[, 2])
env_table[, 1]  <- as.numeric(env_table[, 1])

### env_table[, 1:8]  <- as.numeric(env_table[, 1:8])

#Transforming to log scale 

env_table_scale = scale(log(env_table), center=TRUE,scale = TRUE) # no missing data allowed
env_table_scale
env_table_scale<- data.frame(env_table_scale) #must be transformed to data frame again 

data=data.frame(env_table_scale, meta_table$Site, meta_table$Type, meta_table$Host, row.names=rownames(meta_table))
data

colnames(data) [4] <- "Site"
colnames(data) [5] <- "Type"
colnames(data) [6] <- "Host"
data

# Distance-based redundancy analysis (dbRDA) is an ordination method similar to Redundancy Analysis (rda),
# but it allows non-Euclidean dissimilarity indices, such as Manhattan or Bray-Curtis, Raup-Crick distance. 
# Despite this non-Euclidean feature, the analysis is strictly linear and metric. 

binary_table = t (binary_table)
colnames(binary_table)
raup <- vegdist(binary_table,method="raup")
head(raup)

# null model
raup.min<-capscale (raup ~ 1, data)
plot(raup.min)


# model including all env variables - dbrda(raup ~ ., data) is the same!
raup.all <-dbrda (raup ~ Host + Site +Type, data)
plot(raup.all)

# OrdiR2step: To select the significant explanatory variables
# The criteria for including a variable in the model is based on adjusted variation (R2adj) and Pvalues
# explained by the selected variables and their comparison with R2adj explained by the global model (with all variables);
# if the new variable is not significant or the R2adj of the model including this new variable would exceed the R2adj of the global model,
# the selection will be stopped.

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction =  "both") 
RsquareAdj (raup.min)$adj.r.squared # adjusted R2 explained by all variables


raup.db # new model with selected variables based on on adjusted variation: 
plot(raup.db)
summary(raup.db)
RsquareAdj (raup.db)$adj.r.squared
anova(raup.db)
anova(raup.db, by="term") 
anova(raup.db, by="axis")


mod <- varpart(raup, ~Site, ~Type, ~Host, data = data)
mod
showvarparts(3)
plot(mod, cutoff = -Inf, Xnames = c('Site', 'Type', 'Host'))




# Only no mycorrhizal

nomyc <- subset_taxa(texcoco.binary.models, Myc %in% "t__nomyc") 
nomyc

binary_table <- otu_table(nomyc) #just otu table
binary_table

meta_table <- sample_data(nomyc) #only sample data 
meta_table

Host <- c("Host") 
Site <- c("Site")
Type <- c("Type")
HST <- c("Host", "Site", "Type")


#Transform columns to numeric values, each column separately 

env_table <- subset(meta_table, select = HST)
env_table <-data.frame(env_table) #use as data frame to be able to transform to numeric
env_table[, 3]  <- as.numeric(env_table[, 3])
env_table[, 2]  <- as.numeric(env_table[, 2])
env_table[, 1]  <- as.numeric(env_table[, 1])


#Transforming to log scale 

env_table_scale = scale(log(env_table), center=TRUE,scale = TRUE) # no missing data allowed
env_table_scale
env_table_scale<- data.frame(env_table_scale) #must be transformed to data frame again 

data=data.frame(env_table_scale, meta_table$Site, meta_table$Type, meta_table$Host, row.names=rownames(meta_table))
data

colnames(data) [4] <- "Site"
colnames(data) [5] <- "Type"
colnames(data) [6] <- "Host"
data

# Distance-based redundancy analysis (dbRDA) is an ordination method similar to Redundancy Analysis (rda),
# but it allows non-Euclidean dissimilarity indices, such as Manhattan or Bray-Curtis, Raup-Crick distance. 
# Despite this non-Euclidean feature, the analysis is strictly linear and metric. 

binary_table = t (binary_table)
colnames(binary_table)
raup <- vegdist(binary_table,method="raup")
head(raup)

# null model
raup.min<-capscale (raup ~ 1, data)
plot(raup.min)


# model including all env variables - dbrda(raup ~ ., data) is the same!
raup.all <-dbrda (raup ~ Host + Site +Type, data)
plot(raup.all)



# OrdiR2step: To select the significant explanatory variables
# The criteria for including a variable in the model is based on adjusted variation (R2adj) and Pvalues
# explained by the selected variables and their comparison with R2adj explained by the global model (with all variables);
# if the new variable is not significant or the R2adj of the model including this new variable would exceed the R2adj of the global model,
# the selection will be stopped.

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction =  "both") 
RsquareAdj (raup.min)$adj.r.squared # adjusted R2 explained by all variables


raup.db # new model with selected variables based on on adjusted variation: 
plot(raup.db)
summary(raup.db)
RsquareAdj (raup.db)$adj.r.squared
anova(raup.db)
anova(raup.db, by="term") 
anova(raup.db, by="axis")


mod <- varpart(raup, ~Site, ~Host, ~Type, data = data)
mod
showvarparts(3)
plot(mod, cutoff = -Inf, Xnames = c('Site', 'Type', 'Host'))




# Only ectomycorrhizal samples "ecm" 
ecm <- subset_taxa(texcoco.binary.models, Trophic %in% "a__ecm") 
ecm

binary_table <- otu_table(ecm) #just otu table
binary_table

meta_table <- sample_data(ecm) #only sample data 
meta_table

Host <- c("Host") 
Site <- c("Site")
Type <- c("Type")
HST <- c("Host", "Site", "Type")


#Transform columns to numeric values, each column separately 

env_table <- subset(meta_table, select = HST)
env_table <-data.frame(env_table) #use as data frame to be able to transform to numeric
env_table[, 3]  <- as.numeric(env_table[, 3])
env_table[, 2]  <- as.numeric(env_table[, 2])
env_table[, 1]  <- as.numeric(env_table[, 1])

#Transforming to log scale 

env_table_scale = scale(log(env_table), center=TRUE,scale = TRUE) # no missing data allowed
env_table_scale
env_table_scale<- data.frame(env_table_scale) #must be transformed to data frame again 

data=data.frame(env_table_scale, meta_table$Site, meta_table$Type, meta_table$Host, row.names=rownames(meta_table))
data

colnames(data) [4] <- "Site"
colnames(data) [5] <- "Type"
colnames(data) [6] <- "Host"
data

# Distance-based redundancy analysis (dbRDA) is an ordination method similar to Redundancy Analysis (rda),
# but it allows non-Euclidean dissimilarity indices, such as Manhattan or Bray-Curtis, Raup-Crick distance. 
# Despite this non-Euclidean feature, the analysis is strictly linear and metric. 

binary_table = t (binary_table)
colnames(binary_table)
raup <- vegdist(binary_table,method="raup")
head(raup)

# null model
raup.min<-capscale (raup ~ 1, data)
plot(raup.min)


# model including all env variables - dbrda(raup ~ ., data) is the same!
raup.all <-dbrda (raup ~ Host + Site + Type, data)
plot(raup.all)

# OrdiR2step: To select the significant explanatory variables
# The criteria for including a variable in the model is based on adjusted variation (R2adj) and Pvalues
# explained by the selected variables and their comparison with R2adj explained by the global model (with all variables);
# if the new variable is not significant or the R2adj of the model including this new variable would exceed the R2adj of the global model,
# the selection will be stopped.

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction =  "both") 
RsquareAdj (raup.min)$adj.r.squared # adjusted R2 explained by all variables


raup.db # new model with selected variables based on on adjusted variation: 
plot(raup.db)
summary(raup.db)
RsquareAdj (raup.db)$adj.r.squared
anova(raup.db)
anova(raup.db, by="term") 
anova(raup.db, by="axis")



#Only AM samples "am" 
am <- subset_taxa(texcoco.binary.models, Trophic %in% "a__am") 
am

binary_table <- otu_table(am) #just otu table
binary_table

meta_table <- sample_data(am) #only sample data 
meta_table

Host <- c("Host") 
Site <- c("Site")
Type <- c("Type")
HST <- c("Host", "Site", "Type")


#Transform columns to numeric values, each column separately 

env_table <- subset(meta_table, select = HST)
env_table <-data.frame(env_table) #use as data frame to be able to transform to numeric
env_table[, 3]  <- as.numeric(env_table[, 3])
env_table[, 2]  <- as.numeric(env_table[, 2])
env_table[, 1]  <- as.numeric(env_table[, 1])

#Transforming to log scale 

env_table_scale = scale(log(env_table), center=TRUE,scale = TRUE) # no missing data allowed
env_table_scale
env_table_scale<- data.frame(env_table_scale) #must be transformed to data frame again 

data=data.frame(env_table_scale, meta_table$Site, meta_table$Type, meta_table$Host, row.names=rownames(meta_table))
data

colnames(data) [4] <- "Site"
colnames(data) [5] <- "Type"
colnames(data) [6] <- "Host"
data

# Distance-based redundancy analysis (dbRDA) is an ordination method similar to Redundancy Analysis (rda),
# but it allows non-Euclidean dissimilarity indices, such as Manhattan or Bray-Curtis, Raup-Crick distance. 
# Despite this non-Euclidean feature, the analysis is strictly linear and metric. 

binary_table = t (binary_table)
colnames(binary_table)
raup <- vegdist(binary_table,method="raup")
head(raup)

# null model
raup.min<-capscale (raup ~ 1, data)
plot(raup.min)


# model including all env variables - dbrda(raup ~ ., data) is the same!
raup.all <-dbrda (raup ~ Host + Site + Type, data)
plot(raup.all)

# OrdiR2step: To select the significant explanatory variables
# The criteria for including a variable in the model is based on adjusted variation (R2adj) and Pvalues
# explained by the selected variables and their comparison with R2adj explained by the global model (with all variables);
# if the new variable is not significant or the R2adj of the model including this new variable would exceed the R2adj of the global model,
# the selection will be stopped.

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction =  "both") 
RsquareAdj (raup.min)$adj.r.squared # adjusted R2 explained by all variables


raup.db # new model with selected variables based on on adjusted variation: 
plot(raup.db)
summary(raup.db)
RsquareAdj (raup.db)$adj.r.squared
anova(raup.db)
anova(raup.db, by="term") 
anova(raup.db, by="axis")


mod <- varpart(raup, ~Site, ~Host, ~Type, data = data)
mod
showvarparts(3)
plot(mod, cutoff = -Inf, Xnames = c('Site', 'Type', 'Host'))





#### Models: Texcoco trophic modes in Quercus #### 

#Only no mycorrhizal in Quercus

nomyc <- subset_taxa(texcoco.binary.models, Myc %in% "t__nomyc")
nomyc
nomycq <- subset_samples(nomyc, Host %in% "Quercus")
nomycq

binary_table <- otu_table(nomycq) #just otu table
binary_table

meta_table <- sample_data(nomycq) #only sample data 
meta_table

Host <- c("Host") 
Site <- c("Site")
Type <- c("Type")
ST <- c("Site", "Type")


#Transform columns to numeric values, each column separately 

env_table <- subset(meta_table, select = ST)
env_table <-data.frame(env_table) #use as data frame to be able to transform to numeric
env_table[, 2]  <- as.numeric(env_table[, 2])
env_table[, 1]  <- as.numeric(env_table[, 1])

#Transforming to log scale 

env_table_scale = scale(log(env_table), center=TRUE,scale = TRUE) # no missing data allowed
env_table_scale
env_table_scale<- data.frame(env_table_scale) #must be transformed to data frame again 

data=data.frame(env_table_scale, meta_table$Site, meta_table$Type, row.names=rownames(meta_table))
data

colnames(data) [3] <- "Site"
colnames(data) [4] <- "Type"
data

# Distance-based redundancy analysis (dbRDA) is an ordination method similar to Redundancy Analysis (rda),
# but it allows non-Euclidean dissimilarity indices, such as Manhattan or Bray-Curtis, Raup-Crick distance. 
# Despite this non-Euclidean feature, the analysis is strictly linear and metric. 

binary_table = t (binary_table)
colnames(binary_table)
raup <- vegdist(binary_table,method="raup")
head(raup)

# null model
raup.min<-capscale (raup ~ 1, data)
plot(raup.min)


# model including all env variables - dbrda(raup ~ ., data) is the same!
raup.all <-dbrda (raup ~ Site + Type, data)
plot(raup.all)

# OrdiR2step: To select the significant explanatory variables
# The criteria for including a variable in the model is based on adjusted variation (R2adj) and Pvalues
# explained by the selected variables and their comparison with R2adj explained by the global model (with all variables);
# if the new variable is not significant or the R2adj of the model including this new variable would exceed the R2adj of the global model,
# the selection will be stopped.

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction =  "both") 
RsquareAdj (raup.min)$adj.r.squared # adjusted R2 explained by all variables


raup.db # new model with selected variables based on on adjusted variation: 
plot(raup.db)
summary(raup.db)
RsquareAdj (raup.db)$adj.r.squared
anova(raup.db)
anova(raup.db, by="term") 
anova(raup.db, by="axis")

mod <- varpart(raup, ~Site, ~Type, data = data)
mod
showvarparts(2)
plot(mod, cutoff = -Inf, Xnames = c('Site', 'Type'))




#Only ectomycorrhizal samples "ecm" in Quercus
ecm <- subset_taxa(texcoco.binary.models, Trophic %in% "a__ecm")
ecmq <- subset_samples(ecm, Host %in% "Quercus")
ecmq

binary_table <- otu_table(ecmq) #just otu table
binary_table

meta_table <- sample_data(ecmq) #only sample data 
meta_table

Host <- c("Host") 
Site <- c("Site")
Type <- c("Type")
ST <- c("Site", "Type")


#Transform columns to numeric values, each column separately 

env_table <- subset(meta_table, select = ST)
env_table <-data.frame(env_table) #use as data frame to be able to transform to numeric
env_table[, 2]  <- as.numeric(env_table[, 2])
env_table[, 1]  <- as.numeric(env_table[, 1])

#Transforming to log scale 

env_table_scale = scale(log(env_table), center=TRUE,scale = TRUE) # no missing data allowed
env_table_scale
env_table_scale<- data.frame(env_table_scale) #must be transformed to data frame again 

data=data.frame(env_table_scale, meta_table$Site, meta_table$Type, row.names=rownames(meta_table))
data

colnames(data) [3] <- "Site"
colnames(data) [4] <- "Type"
data

# Distance-based redundancy analysis (dbRDA) is an ordination method similar to Redundancy Analysis (rda),
# but it allows non-Euclidean dissimilarity indices, such as Manhattan or Bray-Curtis, Raup-Crick distance. 
# Despite this non-Euclidean feature, the analysis is strictly linear and metric. 

binary_table = t (binary_table)
colnames(binary_table)
raup <- vegdist(binary_table,method="raup")
head(raup)

# null model
raup.min<-capscale (raup ~ 1, data)
plot(raup.min)


# model including all env variables - dbrda(raup ~ ., data) is the same!
raup.all <-dbrda (raup ~ Site + Type, data)
plot(raup.all)

# OrdiR2step: To select the significant explanatory variables
# The criteria for including a variable in the model is based on adjusted variation (R2adj) and Pvalues
# explained by the selected variables and their comparison with R2adj explained by the global model (with all variables);
# if the new variable is not significant or the R2adj of the model including this new variable would exceed the R2adj of the global model,
# the selection will be stopped.

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction =  "both") 
RsquareAdj (raup.min)$adj.r.squared # adjusted R2 explained by all variables


raup.db # new model with selected variables based on on adjusted variation: 
plot(raup.db)
summary(raup.db)
RsquareAdj (raup.db)$adj.r.squared
anova(raup.db)
anova(raup.db, by="term") 
anova(raup.db, by="axis")

mod <- varpart(raup, ~Site, ~Type, data = data)
mod
showvarparts(2)
plot(mod, cutoff = -Inf, Xnames = c('Site', 'Type'))



#Only AM samples "am" in Quercus
am <- subset_taxa(texcoco.binary.models, Trophic %in% "a__am")
am
amq <- subset_samples(am, Host %in% "Quercus")
amq

binary_table <- otu_table(amq) #just otu table
binary_table

meta_table <- sample_data(amq) #only sample data 
meta_table

Host <- c("Host") 
Site <- c("Site")
Type <- c("Type")
ST <- c("Site", "Type")


#Transform columns to numeric values, each column separately 

env_table <- subset(meta_table, select = ST)
env_table <-data.frame(env_table) #use as data frame to be able to transform to numeric
env_table[, 2]  <- as.numeric(env_table[, 2])
env_table[, 1]  <- as.numeric(env_table[, 1])

#Transforming to log scale 

env_table_scale = scale(log(env_table), center=TRUE,scale = TRUE) # no missing data allowed
env_table_scale
env_table_scale<- data.frame(env_table_scale) #must be transformed to data frame again 

data=data.frame(env_table_scale, meta_table$Site, meta_table$Type, row.names=rownames(meta_table))
data

colnames(data) [3] <- "Site"
colnames(data) [4] <- "Type"
data

# Distance-based redundancy analysis (dbRDA) is an ordination method similar to Redundancy Analysis (rda),
# but it allows non-Euclidean dissimilarity indices, such as Manhattan or Bray-Curtis, Raup-Crick distance. 
# Despite this non-Euclidean feature, the analysis is strictly linear and metric. 

binary_table = t (binary_table)
colnames(binary_table)
raup <- vegdist(binary_table,method="raup")
head(raup)

# null model
raup.min<-capscale (raup ~ 1, data)
plot(raup.min)


# model including all env variables - dbrda(raup ~ ., data) is the same!
raup.all <-dbrda (raup ~ Site + Type, data)
plot(raup.all)

# OrdiR2step: To select the significant explanatory variables
# The criteria for including a variable in the model is based on adjusted variation (R2adj) and Pvalues
# explained by the selected variables and their comparison with R2adj explained by the global model (with all variables);
# if the new variable is not significant or the R2adj of the model including this new variable would exceed the R2adj of the global model,
# the selection will be stopped.

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction =  "both") 
RsquareAdj (raup.min)$adj.r.squared # adjusted R2 explained by all variables


raup.db # new model with selected variables based on on adjusted variation: 
plot(raup.db)
summary(raup.db)
RsquareAdj (raup.db)$adj.r.squared
anova(raup.db)
anova(raup.db, by="term") 
anova(raup.db, by="axis")

mod <- varpart(raup, ~Site, ~Type, data = data)
mod
showvarparts(2)
plot(mod, cutoff = -Inf, Xnames = c('Site', 'Type'))



# Only saprobes in quercus

sap <- subset_taxa(texcoco.binary.models, Trophic %in% "a__sap")
sap
sapq <- subset_samples(sap, Host %in% "Quercus")
sapq

binary_table <- otu_table(sapq) #just otu table
binary_table

meta_table <- sample_data(sapq) #only sample data 
meta_table

Host <- c("Host") 
Site <- c("Site")
Type <- c("Type")
ST <- c("Site", "Type")
HST <- c("Site", "Host", "Type")


#Transform columns to numeric values, each column separately 

env_table <- subset(meta_table, select = ST)
env_table <-data.frame(env_table) #use as data frame to be able to transform to numeric
env_table[, 2]  <- as.numeric(env_table[, 2])
env_table[, 1]  <- as.numeric(env_table[, 1])

#Transforming to log scale 

env_table_scale = scale(log(env_table), center=TRUE,scale = TRUE) # no missing data allowed
env_table_scale
env_table_scale<- data.frame(env_table_scale) #must be transformed to data frame again 

data=data.frame(env_table_scale, meta_table$Site, meta_table$Type, row.names=rownames(meta_table))
data

colnames(data) [3] <- "Site"
colnames(data) [4] <- "Type"
data

# Distance-based redundancy analysis (dbRDA) is an ordination method similar to Redundancy Analysis (rda),
# but it allows non-Euclidean dissimilarity indices, such as Manhattan or Bray-Curtis, Raup-Crick distance. 
# Despite this non-Euclidean feature, the analysis is strictly linear and metric. 

binary_table = t (binary_table)
colnames(binary_table)
raup <- vegdist(binary_table,method="raup")
head(raup)

# null model
raup.min<-capscale (raup ~ 1, data)
plot(raup.min)


# model including all env variables - dbrda(raup ~ ., data) is the same!
raup.all <-dbrda (raup ~ Site + Type, data)
plot(raup.all)

# OrdiR2step: To select the significant explanatory variables
# The criteria for including a variable in the model is based on adjusted variation (R2adj) and Pvalues
# explained by the selected variables and their comparison with R2adj explained by the global model (with all variables);
# if the new variable is not significant or the R2adj of the model including this new variable would exceed the R2adj of the global model,
# the selection will be stopped.

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction =  "both") 
RsquareAdj (raup.min)$adj.r.squared # adjusted R2 explained by all variables


raup.db # new model with selected variables based on on adjusted variation: 
plot(raup.db)
summary(raup.db)
RsquareAdj (raup.db)$adj.r.squared
anova(raup.db)
anova(raup.db, by="term") 
anova(raup.db, by="axis")

mod <- varpart(raup, ~Site, ~Type, data = data)
mod
showvarparts(2)
plot(mod, cutoff = -Inf, Xnames = c('Site', 'Type'))

































#### Models: Texcoco trophic modes in Juniperus ####


# Only nomyc in Juniperus

nomycj <- subset_samples(nomyc, Host %in% "Juniperus")
nomycj

binary_table <- otu_table(nomycj) #just otu table
binary_table

meta_table <- sample_data(nomycj) #only sample data 
meta_table

Host <- c("Host") 
Site <- c("Site")
Type <- c("Type")
ST <- c("Site", "Type")


#Transform columns to numeric values, each column separately 

env_table <- subset(meta_table, select = ST)
env_table <-data.frame(env_table) #use as data frame to be able to transform to numeric
env_table[, 2]  <- as.numeric(env_table[, 2])
env_table[, 1]  <- as.numeric(env_table[, 1])

#Transforming to log scale 

env_table_scale = scale(log(env_table), center=TRUE,scale = TRUE) # no missing data allowed
env_table_scale
env_table_scale<- data.frame(env_table_scale) #must be transformed to data frame again 

data=data.frame(env_table_scale, meta_table$Site, meta_table$Type, row.names=rownames(meta_table))
data

colnames(data) [3] <- "Site"
colnames(data) [4] <- "Type"
data

# Distance-based redundancy analysis (dbRDA) is an ordination method similar to Redundancy Analysis (rda),
# but it allows non-Euclidean dissimilarity indices, such as Manhattan or Bray-Curtis, Raup-Crick distance. 
# Despite this non-Euclidean feature, the analysis is strictly linear and metric. 

binary_table = t (binary_table)
colnames(binary_table)
raup <- vegdist(binary_table,method="raup")
head(raup)

# null model
raup.min<-capscale (raup ~ 1, data)
plot(raup.min)


# model including all env variables - dbrda(raup ~ ., data) is the same!
raup.all <-dbrda (raup ~ Site + Type, data)
plot(raup.all)

# OrdiR2step: To select the significant explanatory variables
# The criteria for including a variable in the model is based on adjusted variation (R2adj) and Pvalues
# explained by the selected variables and their comparison with R2adj explained by the global model (with all variables);
# if the new variable is not significant or the R2adj of the model including this new variable would exceed the R2adj of the global model,
# the selection will be stopped.

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction =  "both") 
RsquareAdj (raup.min)$adj.r.squared # adjusted R2 explained by all variables


raup.db # new model with selected variables based on on adjusted variation: 
plot(raup.db)
summary(raup.db)
RsquareAdj (raup.db)$adj.r.squared
anova(raup.db)
anova(raup.db, by="term") 
anova(raup.db, by="axis")

mod <- varpart(raup, ~Site, ~Type, data = data)
mod
showvarparts(2)
plot(mod, cutoff = -Inf, Xnames = c('Site', 'Type'))



#Only am samples in Juniperus 
am <- subset_taxa(texcoco.binary.models, Trophic %in% "a__am")
amj <- subset_samples(am, Host %in% "Juniperus")
amj

binary_table <- otu_table(amj) #just otu table
binary_table

meta_table <- sample_data(amj) #only sample data 
meta_table

Host <- c("Host") 
Site <- c("Site")
Type <- c("Type")
ST <- c("Site", "Type")


#Transform columns to numeric values, each column separately 

env_table <- subset(meta_table, select = ST)
env_table <-data.frame(env_table) #use as data frame to be able to transform to numeric
env_table[, 2]  <- as.numeric(env_table[, 2])
env_table[, 1]  <- as.numeric(env_table[, 1])

#Transforming to log scale 

env_table_scale = scale(log(env_table), center=TRUE,scale = TRUE) # no missing data allowed
env_table_scale
env_table_scale<- data.frame(env_table_scale) #must be transformed to data frame again 

data=data.frame(env_table_scale, meta_table$Site, meta_table$Type, row.names=rownames(meta_table))
data

colnames(data) [3] <- "Site"
colnames(data) [4] <- "Type"
data

# Distance-based redundancy analysis (dbRDA) is an ordination method similar to Redundancy Analysis (rda),
# but it allows non-Euclidean dissimilarity indices, such as Manhattan or Bray-Curtis, Raup-Crick distance. 
# Despite this non-Euclidean feature, the analysis is strictly linear and metric. 

binary_table = t (binary_table)
colnames(binary_table)
raup <- vegdist(binary_table,method="raup")
head(raup)

# null model
raup.min<-capscale (raup ~ 1, data)
plot(raup.min)


# model including all env variables - dbrda(raup ~ ., data) is the same!
raup.all <-dbrda (raup ~ Site + Type, data)
plot(raup.all)

# OrdiR2step: To select the significant explanatory variables
# The criteria for including a variable in the model is based on adjusted variation (R2adj) and Pvalues
# explained by the selected variables and their comparison with R2adj explained by the global model (with all variables);
# if the new variable is not significant or the R2adj of the model including this new variable would exceed the R2adj of the global model,
# the selection will be stopped.

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction =  "both") 
RsquareAdj (raup.min)$adj.r.squared # adjusted R2 explained by all variables


raup.db # new model with selected variables based on on adjusted variation: 
plot(raup.db)
summary(raup.db)
RsquareAdj (raup.db)$adj.r.squared
anova(raup.db)
anova(raup.db, by="term") 
anova(raup.db, by="axis")

mod <- varpart(raup, ~Site, ~Type, data = data)
mod
showvarparts(2)
plot(mod, cutoff = -Inf, Xnames = c('Site', 'Type'))


# Only ecm in juniperus

ecm
ecmj <- subset_samples(ecm, Host %in% "Juniperus")
ecmj

binary_table <- otu_table(ecmj) #just otu table
binary_table

meta_table <- sample_data(ecmj) #only sample data 
meta_table

Host <- c("Host") 
Site <- c("Site")
Type <- c("Type")
ST <- c("Site", "Type")


#Transform columns to numeric values, each column separately 

env_table <- subset(meta_table, select = ST)
env_table <-data.frame(env_table) #use as data frame to be able to transform to numeric
env_table[, 2]  <- as.numeric(env_table[, 2])
env_table[, 1]  <- as.numeric(env_table[, 1])

#Transforming to log scale 

env_table_scale = scale(log(env_table), center=TRUE,scale = TRUE) # no missing data allowed
env_table_scale
env_table_scale<- data.frame(env_table_scale) #must be transformed to data frame again 

data=data.frame(env_table_scale, meta_table$Site, meta_table$Type, row.names=rownames(meta_table))
data

colnames(data) [3] <- "Site"
colnames(data) [4] <- "Type"
data

# Distance-based redundancy analysis (dbRDA) is an ordination method similar to Redundancy Analysis (rda),
# but it allows non-Euclidean dissimilarity indices, such as Manhattan or Bray-Curtis, Raup-Crick distance. 
# Despite this non-Euclidean feature, the analysis is strictly linear and metric. 

binary_table = t (binary_table)
colnames(binary_table)
raup <- vegdist(binary_table,method="raup")
head(raup)

# null model
raup.min<-capscale (raup ~ 1, data)
plot(raup.min)


# model including all env variables - dbrda(raup ~ ., data) is the same!
raup.all <-dbrda (raup ~ Site + Type, data)
plot(raup.all)

# OrdiR2step: To select the significant explanatory variables
# The criteria for including a variable in the model is based on adjusted variation (R2adj) and Pvalues
# explained by the selected variables and their comparison with R2adj explained by the global model (with all variables);
# if the new variable is not significant or the R2adj of the model including this new variable would exceed the R2adj of the global model,
# the selection will be stopped.

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction =  "both") 
RsquareAdj (raup.min)$adj.r.squared # adjusted R2 explained by all variables


raup.db # new model with selected variables based on on adjusted variation: 
plot(raup.db)
summary(raup.db)
RsquareAdj (raup.db)$adj.r.squared
anova(raup.db)
anova(raup.db, by="term") 
anova(raup.db, by="axis")

mod <- varpart(raup, ~Site, ~Type, data = data)
mod
showvarparts(2)
plot(mod, cutoff = -Inf, Xnames = c('Site', 'Type'))



# Only saprobes in juniperus

sap <- subset_taxa(texcoco.binary.models, Trophic %in% "a__sap")
sap
sapj <- subset_samples(sap, Host %in% "Juniperus")
sapj

binary_table <- otu_table(sapj) #just otu table
binary_table

meta_table <- sample_data(sapj) #only sample data 
meta_table

Host <- c("Host") 
Site <- c("Site")
Type <- c("Type")
ST <- c("Site", "Type")
HST <- c("Site", "Host", "Type")


#Transform columns to numeric values, each column separately 

env_table <- subset(meta_table, select = ST)
env_table <-data.frame(env_table) #use as data frame to be able to transform to numeric
env_table[, 2]  <- as.numeric(env_table[, 2])
env_table[, 1]  <- as.numeric(env_table[, 1])

#Transforming to log scale 

env_table_scale = scale(log(env_table), center=TRUE,scale = TRUE) # no missing data allowed
env_table_scale
env_table_scale<- data.frame(env_table_scale) #must be transformed to data frame again 

data=data.frame(env_table_scale, meta_table$Site, meta_table$Type, row.names=rownames(meta_table))
data

colnames(data) [3] <- "Site"
colnames(data) [4] <- "Type"
data

# Distance-based redundancy analysis (dbRDA) is an ordination method similar to Redundancy Analysis (rda),
# but it allows non-Euclidean dissimilarity indices, such as Manhattan or Bray-Curtis, Raup-Crick distance. 
# Despite this non-Euclidean feature, the analysis is strictly linear and metric. 

binary_table = t (binary_table)
colnames(binary_table)
raup <- vegdist(binary_table,method="raup")
head(raup)

# null model
raup.min<-capscale (raup ~ 1, data)
plot(raup.min)


# model including all env variables - dbrda(raup ~ ., data) is the same!
raup.all <-dbrda (raup ~ Site + Type, data)
plot(raup.all)

# OrdiR2step: To select the significant explanatory variables
# The criteria for including a variable in the model is based on adjusted variation (R2adj) and Pvalues
# explained by the selected variables and their comparison with R2adj explained by the global model (with all variables);
# if the new variable is not significant or the R2adj of the model including this new variable would exceed the R2adj of the global model,
# the selection will be stopped.

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction =  "both") 
RsquareAdj (raup.min)$adj.r.squared # adjusted R2 explained by all variables


raup.db # new model with selected variables based on on adjusted variation: 
plot(raup.db)
summary(raup.db)
RsquareAdj (raup.db)$adj.r.squared
anova(raup.db)
anova(raup.db, by="term") 
anova(raup.db, by="axis")

mod <- varpart(raup, ~Site, ~Type, data = data)
mod
showvarparts(2)
plot(mod, cutoff = -Inf, Xnames = c('Site', 'Type'))








#### Nested models ####

nested.anova.dbrda(raup, data) #try?


nested.npmanova(formula, data, method="euc", permutations=100, warnings=FALSE) #try? 















#### IZTA ####

soil <- subset_samples(izta.binary.models, Type %in% "soil") # phyloseq object subset with only soil samples
soil


binary_table <- otu_table(soil) #just otu table
binary_table

meta_table <- sample_data(soil) #only sample data and only from soil samples
meta_table

env <- c("pH","Pdis","Ca","Mg","K","Na","H","SoilM") # Al removed becuase of missing data
micro <- c("Ca","Mg","K","Na","H")

#Transform columns to numeric values, each column separately 

env_table <- subset(meta_table, select = env)
env_table <-data.frame(env_table) #use as data frame to be able to transform to numeric
env_table[, 8]  <- as.numeric(env_table[, 8])
env_table[, 7]  <- as.numeric(env_table[, 7])
env_table[, 6]  <- as.numeric(env_table[, 6])
env_table[, 5]  <- as.numeric(env_table[, 5])
env_table[, 4]  <- as.numeric(env_table[, 4])
env_table[, 3]  <- as.numeric(env_table[, 3])
env_table[, 2]  <- as.numeric(env_table[, 2])
env_table[, 1]  <- as.numeric(env_table[, 1])

### env_table[, 1:8]  <- as.numeric(env_table[, 1:8])

#Transforming to log scale 

env_table_scale = scale(log(env_table), center=TRUE,scale = TRUE) # no missing data allowed
env_table_scale
env_table_scale<- data.frame(env_table_scale) #must be transformed to data frame again 

data=data.frame(env_table_scale, meta_table$Site, meta_table$Host, row.names=rownames(meta_table))
data
colnames(data) [9] <- "Site"
colnames(data) [10] <- "Host"
data

micro_table <- subset(meta_table, select = micro)
micro_table

# Distance-based redundancy analysis (dbRDA) is an ordination method similar to Redundancy Analysis (rda),
# but it allows non-Euclidean dissimilarity indices, such as Manhattan or Bray-Curtis, Raup-Crick distance. 
# Despite this non-Euclidean feature, the analysis is strictly linear and metric. 

binary_table = t (binary_table)
colnames(binary_table)
raup <- vegdist(binary_table,method="raup")
head(raup)

# null model
raup.min<-capscale (raup ~ 1, env_table_scale)
plot(raup.min)

# model including all env variables - dbrda(raup ~ ., data) is the same!
raup.all <-dbrda (raup ~ pH + Pdis + Ca + Mg + K + Na + H + SoilM, env_table_scale) 

# hypothesis model
raup.test<-dbrda(raup ~ pH + Pdis, env_table_scale)

# which is the best model?
# if both are significant, you choose the simplest model (raup.test)
anova (raup.all, raup.test)


# OrdiR2step: To select the significant explanatory variables
# The criteria for including a variable in the model is based on adjusted variation (R2adj) and Pvalues
# explained by the selected variables and their comparison with R2adj explained by the global model (with all variables);
# if the new variable is not significant or the R2adj of the model including this new variable would exceed the R2adj of the global model,
# the selection will be stopped.

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction =  "both")
RsquareAdj (raup.min)$adj.r.squared # adjusted R2 explained by all variables


raup.db # new model with selected variables based on on adjusted variation: 
plot(raup.db)
summary(raup.db)
RsquareAdj (raup.db)$adj.r.squared
anova(raup.db)
anova(raup.db, by="term") 
anova(raup.db, by="axis")


# model including site
raup.min<-capscale (raup ~ 1, data)
raup.all <-dbrda (raup ~ pH + Pdis + Ca + Mg + K + Na + H + Site, data) 
# removed Pdis, error message: Error in dimnames(u) <- list(dnam[[1]], c(axnam, negnam)) : length of 'dimnames' [2] not equal to array extent

# OrdiR2step: To select the significant explanatory variables
# The criteria for including a variable in the model is based on adjusted variation (R2adj) and Pvalues
# explained by the selected variables and their comparison with R2adj explained by the global model (with all variables);
# if the new variable is not significant or the R2adj of the model including this new variable would exceed the R2adj of the global model,
# the selection will be stopped.

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction =  "both")
RsquareAdj (raup.min)$adj.r.squared # adjusted R2 explained by all variables


raup.db # new model with selected variables based on on adjusted variation:  
plot(raup.db)
summary(raup.db)
RsquareAdj (raup.db)$adj.r.squared
anova(raup.db)
anova(raup.db, by="term") 
anova(raup.db, by="axis")















