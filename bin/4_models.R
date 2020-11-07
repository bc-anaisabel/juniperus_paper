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
library(car)
library(lme4)
library(PerformanceAnalytics)
library(tables)
library(Rmisc)

R.version

# Set working directory to source file


source("../bin/1_Filter_otu_table.R")

juniperus.rel <- phyloseq.rel  #change name of sequence abundance table for this script 

juniperus.rel
# juniperus.rel_OTU2 # TO FIX! (Right now it calculates reads sum and not taxa sum)

#### Subset of binomial tables: Proyecto Texcoco ####

texcoco.binary.models <- subset_samples(binary_table_OTU2, Project %in% "Texcoco")
texcoco.binary.models
sample_data(texcoco.binary.models)

# Remove OTUs that are not present in any of the samples of the subset (OTUs that are only present either in Izta or Texcoco) 
any(taxa_sums(texcoco.binary.models) == 0)
taxa_sums(texcoco.binary.models) [1:10]
texcoco.binary.models <- prune_taxa(taxa_sums(texcoco.binary.models) > 0, texcoco.binary.models)
any(taxa_sums(texcoco.binary.models) == 0)
texcoco.binary.models


#### Models using metadata: Texcoco only soil samples ####

#subset 
soil <- subset_samples(texcoco.binary.models, Type %in% "soil") # phyloseq object subset with only soil samples
soil

binary_table <- otu_table(soil) #just otu table
binary_table

meta_table <- sample_data(soil) #only sample data and only from soil samples
meta_table

env <- c("pH","Pdis","Ca","Mg","K","Na","H","Al","SoilM","CN") 

#Transform columns to numeric values, each column separately 

env_table <- subset(meta_table, select = env)
env_table <-data.frame(env_table) #use as data frame to be able to transform to numeric
env_table[, 10] <- as.numeric(as.character(env_table[, 10]))
env_table[, 9]  <- as.numeric(as.character(env_table[, 9]))
env_table[, 8]  <- as.numeric(as.character(env_table[, 8]))
env_table[, 7]  <- as.numeric(as.character(env_table[, 7]))
env_table[, 6]  <- as.numeric(as.character(env_table[, 6]))
env_table[, 5]  <- as.numeric(as.character(env_table[, 5]))
env_table[, 4]  <- as.numeric(as.character(env_table[, 4]))
env_table[, 3]  <- as.numeric(as.character(env_table[, 3]))
env_table[, 2]  <- as.numeric(as.character(env_table[, 2]))
env_table[, 1]  <- as.numeric(as.character(env_table[, 1]))

### env_table[, 1:8]  <- as.numeric(env_table[, 1:8])

#Transforming to log scale 

env_table_scale = scale(log(env_table), center=TRUE,scale = TRUE) # no missing data allowed
env_table_scale
env_table_scale<- data.frame(env_table_scale) #must be transformed to data frame again 

data=data.frame(env_table_scale, meta_table$Site, meta_table$Host, row.names=rownames(meta_table))
data
colnames(data) [12] <- "Site"
colnames(data) [13] <- "Host"
data

# Create soil table with average and SD values 

data=data.frame(env_table, meta_table$Site, meta_table$Host, row.names=rownames(meta_table))
data
colnames(data) [11] <- "Site"
colnames(data) [12] <- "Host"
data

###summary(data)
# using exported soilvariables data 
data <-read.csv("soilvariables.csv", row.names = 1)

x <- colnames(data[,1:10])

# loop that calculates summary statistics for each variable
# changes the third column name to "mean" instead of variable name
# and adds additional column titled 'Variable' with variable name

for (i in x){
  a <- summarySE(data, measurevar= i, groupvars=c("Site"), na.rm = T)
  names(a)[names(a) == i] <- "mean"
  a$Variable <- i
  assign(i,a)}

Bind_soilvar<-rbind(pH, Pdis, Ca, Mg, K, Na, H, Al, SoilM, CN)
write.csv(Bind_soilvar, col.names=TRUE, file = "bind_soilvar.csv")

# Add results of ANOVA and pairwise Tukey tests in table

aov<-aov (data$CN ~ Site, data = data)
summary(aov)
TukeyHSD(aov)


#### the longer option #### 

pH <- summarySE(data, measurevar= "pH", groupvars=c("Site", "Host"), na.rm = TRUE)
Pdis <-summarySE(data, measurevar= "Pdis", groupvars=c("Site", "Host"), na.rm = TRUE)
Ca <-summarySE(data, measurevar= "Ca", groupvars=c("Site", "Host"), na.rm = TRUE)
Mg <-summarySE(data, measurevar= "Mg", groupvars=c("Site", "Host"), na.rm = TRUE)
K <-summarySE(data, measurevar= "K", groupvars=c("Site", "Host"), na.rm = TRUE)
Na <-summarySE(data, measurevar= "Na", groupvars=c("Site", "Host"), na.rm = TRUE)
H <-summarySE(data, measurevar= "H", groupvars=c("Site", "Host"), na.rm = TRUE)
Al <-summarySE(data, measurevar= "Al", groupvars=c("Site", "Host"), na.rm = TRUE)
SoilM <-summarySE(data, measurevar= "SoilM", groupvars=c("Site", "Host"), na.rm = TRUE)
CN <-summarySE(data, measurevar= "CN", groupvars=c("Site", "Host"), na.rm = TRUE)

soilvar= (data.frame(pH,Pdis,Ca,Mg,K,Na,H,Al,SoilM,CN))
write.csv.tabular(soilvar, file = "soilvar.csv")





#### Scatterplots ####

## Correlation matrix with p-values. See http://goo.gl/nahmV for documentation of this function
cor.prob <- function (X, dfr = nrow(X) - 2) {
  R <- cor(X, use="pairwise.complete.obs")
  above <- row(R) < col(R)
  r2 <- R[above]^2
  Fstat <- r2 * dfr/(1 - r2)
  R[above] <- 1 - pf(Fstat, 1, dfr)
  R[row(R) == col(R)] <- NA
  R
}

## Use this to dump the cor.prob output to a 4 column matrix
## with row/column indices, correlation, and p-value.
## See StackOverflow question: http://goo.gl/fCUcQ
flattenSquareMatrix <- function(m) {
  if( (class(m) != "matrix") | (nrow(m) != ncol(m))) stop("Must be a square matrix.") 
  if(!identical(rownames(m), colnames(m))) stop("Row and column names must be equal.")
  ut <- upper.tri(m)
  data.frame(i = rownames(m)[row(m)[ut]],
             j = rownames(m)[col(m)[ut]],
             cor=t(m)[ut],
             p=m[ut])
}

# get some data from the mtcars built-in dataset
mydata <- data[, c(1:8)]

# correlation matrix
cor(mydata)

# correlation matrix with p-values
cor.prob(mydata)

# "flatten" that table
flattenSquareMatrix(cor.prob(mydata))

# plot the data

chart.Correlation(mydata)


#### Models using site, host and type as factors ####

# subset for different trophic modes and subset by sample type and host 

model<- texcoco.binary.models #from binary_table_OTU2

sample_data(model)

foram = subset_samples(texcoco.binary.models, phinchID != "F3-J5s-2-2621810") #take out samples with no fungi (either ecm or am)

foram

model <- foram

model


any(taxa_sums(model) == 0)

model

#Create tables for modeling using raup 

binary_table <- otu_table(model) 
binary_table

meta_table <- sample_data(model) 
meta_table

Host <- c("Host") 
Site <- c("Site")
Type <- c("Type")
HS <- c("Site", "Host")
ST <- c("Site", "Type")
HST <- c("Site", "Host", "Type")

#Transform columns to numeric values, each column separately 

env_table <- subset(meta_table, select = HST)
env_table <-data.frame(env_table) #use as data frame to be able to transform to numeric
env_table[, 3]  <- as.numeric(env_table[, 3]) #not use this line if subsetting using host (or type)
env_table[, 2]  <- as.numeric(env_table[, 2])
env_table[, 1]  <- as.numeric(env_table[, 1])


#Transforming to log scale 

env_table_scale = scale(log(env_table), center=TRUE,scale = TRUE) # no missing data allowed
env_table_scale
env_table_scale<- data.frame(env_table_scale) #must be transformed to data frame again 

data=data.frame(env_table_scale, meta_table$Site, meta_table$Host, meta_table$Type, row.names=rownames(meta_table)) #add or remove column name as needed 
data
colnames(data) [4] <- "Site"
colnames(data) [5] <- "Host"
colnames(data) [6] <- "Type"
data

binary_table = t (binary_table)
colnames(binary_table)
raup <- vegdist(binary_table,method="raup")
head(raup)
raup
scores(raup)

# null model
raup.min<-capscale (raup ~ 1, data)
plot(raup.min)

# model including host, site and type 
raup.all <-dbrda (raup ~ Host + Site + Type, data)
plot(raup.all)

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction =  "both") 
RsquareAdj (raup.min)$adj.r.squared # adjusted R2 explained by all variables


raup.db # new model with selected variables based on on adjusted variation:  
plot(raup.db)
summary( raup.db)
RsquareAdj (raup.db)$adj.r.squared
anova(raup.db)
anova(raup.db, by="term") 
anova(raup.db, by="axis")

# To construct better plot using ggplot: 



#import to table or matrix? and then to ggplot 

#example?: 
# sample_data(texcoco.binary.models)$Site = factor(sample_data(texcoco.binary.models)$Site, levels=c("native","mixed","perturbated"))
# NMDS Bray Texcoco
#bray_nmds = distance(subset.texcoco.binary.beta, method = "bray")
#ordination = ordinate(subset.texcoco.binary.beta, method = "NMDS", distance = bray_nmds)
#p1 <- plot_ordination(subset.texcoco.binary.beta, ordination, color="Host", shape = "Type", title = "") + theme(aspect.ratio=1)+geom_point(size=3)
#print(p1)
#p1 + facet_wrap(~Site)

# Variation partitioning 

mod <- varpart(raup, data$Site, data$Type, data$Host)
mod

#Use fill colours
showvarparts(3, bg = 2:4)
plot(mod, bg=2:4, Xnames = c('Site', 'Type', 'Host'))


dbrda.result <- dbrda (raup ~ Host + Site + Type, data)
anova(dbrda.result, step=200, perm.max=200)



# Three explanatory matrices 
mod <- varpart(mite, ~ SubsDens + WatrCont, ~ Substrate + Shrub + Topo,
               mite.pcnm, data=mite.env, transfo="hel")
mod
showvarparts(3, bg=2:4)
plot(mod, bg=2:4)
# An alternative formulation of the previous model using
# matrices mm1 amd mm2 and Hellinger transformed species data
mm1 <- model.matrix(~ SubsDens + WatrCont, mite.env)[,-1]
mm2 <- model.matrix(~ Substrate + Shrub + Topo, mite.env)[, -1]
mite.hel <- decostand(mite, "hel")
mod <- varpart(mite.hel, mm1, mm2, mite.pcnm)
# Use RDA to test fraction [a]
# Matrix can be an argument in formula
rda.result <- rda(mite.hel ~ mm1 + Condition(mm2) +
                    Condition(as.matrix(mite.pcnm)))
anova(rda.result, step=200, perm.max=200)










# Diagnostic tools for constrained ordination

vif.cca(raup.db)



# General linear models 

subset.glm <- glm (raup ~ Host, data = data)
summary(subset.glm)













