# Part 7. Script for Illumina MiSeq2x300 Juniperus projects: Texcoco 
# In R this is step 4. Network presence-absence  
# January, 2021

library(devtools)
library(SpiecEasi)
library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(dplyr) # data handling
library(SpiecEasi) # Network analysis for sparse compositional data  
library(network)
library(intergraph)
library(ggnet)
library(igraph)
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


#### Network using sna #### 

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
taxa_names(tax_table(subset))

color<-as.matrix(tax_table(subset))
color<-as.data.frame(color)
color<- color$Family

#network_host_2<-cbind(network_host,color)
#network_host_2

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

#write the sample out of R for issue example

network_df<-write.csv(network_host, file = "network_df.csv")
color<-as.data.frame(color)
#network_host_t<-as.matrix(network_host)
#network_host_tt<-as.network(network_host_t)

#Gplot
gplot(network_host, thresh = 0.2, displaylabels = TRUE, usearrows=FALSE, 
      legend(x=1,y=-1, pch=21, col = "#777777", 
             pt.cex=2, cex=.8, bty="n", ncol=1), vertex.col = color$Family)


gplot(network_host, thresh = 0.2, displaylabels = TRUE, usearrows=FALSE, 
      legend(x=1,y=-1, color, pch=21,  
             pt.cex=2, cex=.8, bty="n", ncol=1), vertex.col = color$Family)


#Other type of plot 

pal2 <-polychrome(27)

par(mfrow=c(1,2), xpd=T)

gplot(as.one.mode(network_host),
      displaylabels = TRUE,
      gmode="graph",
      label.cex=1, vertex.col = color$Family, vertex.cex=1)

palette(polychrome(n=27))
gplot(network_host, gmode="graph", jitter=FALSE,
      label = color$Family, 
      boxed.labels=FALSE, label.pos=1, label.cex=1, vertex.cex=2,
      vertex.col= as.numeric(color$Family))

#Try ggnet 

library(network)
library(sna)
library(ggplot2)
library(GGally)
library(RColorBrewer)
library(intergraph)



net = network(otu_table(subset), directed = FALSE)
ggnet2(net, node.size = 3, edge.size = 1, node.color = 'mode', edge.color = "grey", label = TRUE)

#or

network_host<-t(network_host)
network_host<-as.data.frame(network_host, stringsAsFactors = F)
network_host[, 2] <- as.character(network_host[, 2])
net = network(network_host, directed = FALSE)
ggnet2(net, node.size = 3, edge.size = 1, node.color = "mode", edge.color = "grey", label = TRUE)
read.csv("network_df.csv")
nuevonet<-read.csv("network_df.csv")
library(tidyverse)
nuevonet %>% remove_rownames %>% column_to_rownames(var="Family")




color<-as.character(color)
familias<-levels(color)
color_vector







#"f__ Glomeraceae", "f__ Paraglomeraceae", "f__ Claroideoglomeraceae", "f__", "f__ Acaulosporaceae","f__ Ambisporaceae", "f__ Gigasporaceae"


#### Networks using igraph ####

#Other format
abc<-t(network_host)

#Plot
plotweb(network_host)

#Plot
visweb(network_host)

network_host

######################







head(tax_table(subset))
otu.c <- t(otu_table(subset)@.Data) #extract the otu table from phyloseq object
tax.c <- as.data.frame(tax_table(subset)@.Data)#extract the taxonomy information

head(tax.c)

am.net <- as.network(otu.c)

network::set.edge.attribute(am.net, "color", ifelse(am.net %e% "weight" > 0, "steelblue", "orange"))

colnames(tax_table(subset))

net.c <- spiec.easi(otu.c, method='mb', icov.select.params=list(rep.num=50)) # reps have to increases for real data

class(net.c)

n.c <- symBeta(getOptBeta(net.c))

colnames(n.c) <- rownames(n.c) <- colnames(otu.c)

vsize <- log2(apply(otu.c, 2, mean)) # add log abundance as properties of vertex/nodes.

am.ig <- graph.adjacency(n.c, mode='undirected', add.rownames = TRUE, weighted = TRUE)
am.ig # we can see all the attributes and weights

plot(am.ig)
?layout_with_fr
coords.fdr = layout_with_fr(am.ig)
E(am.ig)[weight > 0]$color<-"steelblue" #now color the edges based on their values positive is steelblue
E(am.ig)[weight < 0]$color<-"orange"  #now color the edges based on their values

plot(am.ig, layout=coords.fdr, vertex.size = 2, vertex.label.cex = 0.5)


am.net <- asNetwork(am.ig)
network::set.edge.attribute(am.net, "color", ifelse(am.net %e% "weight" > 0, "steelblue", "orange"))
colnames(tax_table(subset))
phyla <- map_levels(colnames(otu.c), from = "Species", to = "Myc", tax_table(subset))
am.net %v% "Phylum" <- phyla
am.net %v% "nodesize" <- vsize


phyla <- map_levels(colnames(otu.c), from = "Species", to = "Myc", tax_table(subset))
am.net %v% "Phylum" <- phyla
am.net %v% "nodesize" <- vsize

mycolors <- scale_color_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928"))
p <- ggnet2(am.net, node.color = "Phylum", 
            label = TRUE, node.size = "nodesize", 
            label.size = 2, edge.color = "color") + guides(color=guide_legend(title="Phylum"), size = FALSE) + mycolors

se_mb <- spiec.easi(subset, method='mb', 
                    lambda.min.ratio=1e-2, nlambda=20, 
                    pulsar.params=list(rep.num=50, ncores=3))


e_net <- adj2igraph(getRefit(se_mb), 
                    rmEmptyNodes = TRUE, diag = FALSE, 
                    vertex.attr = list(name = taxa_names(subset))) # Usamos el ID de las taxas para nombrar los vértices o nodos de la red
plot_network(e_net, subset, type = "taxa", color = "Family", shape = "Phylum", label = NULL)

net_class <- as_adjacency_matrix(e_net, type = "both")
net_class <- network(as.matrix(net_class), 
                     vertex.attrnames = sample.names(subset), 
                     matrix.type = "adjacency", directed = F)
ggnet2(net_class)
net <- e_net
nodes <- V(net) # Nodos
edges <- E(net) # Bordes
node.names <- V(net)$name # Nombre de nodos
num.nodes <- vcount(net) # Número total de nodos
num.edges <- ecount(net) # Número total de bordes
# V() es por vértices o nodos
# E() es por edges



ggnet2(net)


imc<-cluster_infomap(network_host)
membership(imc)

plot(g, edge.arrow.size=.2,vertex.label=NA, 
     layout=layout_with_fr,
     vertex.color=imc$membership)
