# Script for Illumina MiSeq2x300 Juniperus project (Texcoco and Izta)
# Part 1. Filtering steps (negative control, relative abundance, 0 presence OTUs, presence/absence)
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


# Import biom data (includes otu table, taxonomy table and sample data)

setwd("~/Documents/Ana+Camille/MiSeq2x300_Dec2019/Juniperus/data/taxonomy")

juniperus <- import_biom("taxonomy.biom")

juniperus

#Data looks like this: 

head(otu_table(juniperus))
head(tax_table(juniperus))
sample_data(juniperus)

# Get colors and theme set 
display.brewer.all(type="seq")
display.brewer.all(type="div")
display.brewer.all(type="qual")

theme_set(theme_bw())

# View and change taxonomic ranks 
rank_names(juniperus) 
colnames(tax_table(juniperus)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species") # change names of taxonomic table
# Not yet: first have to do funguild table : colnames(tax_table(juniperus.guild)) <- c("Guild", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species") # change names of taxonomic table

rank_names(juniperus) # check changes

########################################## Negative control filter

# If necessary, subset data by PCR batch. One negative control for each PCR, so those sequences need to be removed each PCR batch.
PCR_batch = subset_samples(juniperus, PCRbatch == "PCR1") 
PCR_batch

# move PCR_batch to a dataframe to subtract things.
PCR_batch_matrix = (as(otu_table(PCR_batch), "matrix"))
PCR_batch_matrix
PCR_batch_df = as.data.frame(PCR_batch_matrix)
row.names(PCR_batch_df)
colnames(PCR_batch_df)

# Create a vector from the Negative control
NEG1=as.vector(PCR_batch_df$`B7-NEG1-2-2621838`)
NEG1

#Insert the vector into a matrix
test1=cbind(NEG1)
test1

# Replicate the column for the number of samples
test2 = test1[,rep(1,24)]
test2

#Subtract values from original matrix using the matrix based on repetition of negative control 
PCRbatch1_clean<-PCR_batch_matrix-test2
PCRbatch1_clean

# replace all negative numbers with 0, since we can't have negative seqs reads 
PCRbatch1_cleaned <- replace(PCRbatch1_clean, PCRbatch1_clean < 0, 0)
PCRbatch1_cleaned[,1] # check that it worked

# Do those steps again for each PCR batch to remove the contamination from the negative controls. 

#PCR2 Negative 2
PCR_batch2 = subset_samples(juniperus, PCRbatch == "PCR2") 
PCR_batch2

# move PCR_batch to a dataframe to subtract things.
PCR_batch_matrix = (as(otu_table(PCR_batch2), "matrix"))
PCR_batch_matrix
PCR_batch2_df = as.data.frame(PCR_batch_matrix)
row.names(PCR_batch2_df)
colnames(PCR_batch2_df)

# Create a vector from the Negative control
NEG2=as.vector(PCR_batch2_df$`H9-NEG2-2-2621910`)
NEG2

#Insert the vector into a matrix
test1=cbind(NEG2)
test1

# Replicate the column for the number of samples
test2 = test1[,rep(1,20)]
test2

#Subtract values from original matrix using the matrix based on repetition of negative control
PCRbatch2_clean<-PCR_batch_matrix-test2
PCRbatch2_clean

# replace all negative numbers with 0, since we can't have negative seqs reads 
PCRbatch2_cleaned <- replace(PCRbatch2_clean, PCRbatch2_clean < 0, 0)
PCRbatch2_cleaned[,1] # check that it worked

#PCR3 Negative 3

PCR_batch3 = subset_samples(juniperus, PCRbatch == "PCR3") 
PCR_batch3

# move PCR_batch to a dataframe to subtract things.
PCR_batch_matrix = (as(otu_table(PCR_batch3), "matrix"))
PCR_batch_matrix
PCR_batch3_df = as.data.frame(PCR_batch_matrix)
row.names(PCR_batch3_df)
colnames(PCR_batch3_df)

# Create a vector from the Negative control
NEG3=as.vector(PCR_batch3_df$`C2-NEG3-2-2621749`)
NEG3

#Insert the vector into a matrix
test1=cbind(NEG3)
test1

# Replicate the column for the number of samples
test2 = test1[,rep(1,13)]
test2

#Subtract values from original matrix using the matrix based on repetition of negative control
PCRbatch3_clean<-PCR_batch_matrix-test2
PCRbatch3_clean

# replace all negative numbers with 0, since we can't have negative seqs reads 
PCRbatch3_cleaned <- replace(PCRbatch3_clean, PCRbatch3_clean < 0, 0)
PCRbatch3_cleaned[,1] # check that it worked

#PCR4 Negative 4 

PCR_batch4 = subset_samples(juniperus, PCRbatch == "PCR4") 
PCR_batch4

# move PCR_batch to a dataframe to subtract things.
PCR_batch_matrix = (as(otu_table(PCR_batch4), "matrix"))
PCR_batch_matrix
PCR_batch4_df = as.data.frame(PCR_batch_matrix)
row.names(PCR_batch4_df)
colnames(PCR_batch4_df)

# Create a vector from the Negative control
NEG4=as.vector(PCR_batch4_df$`C8-NEG4-2-2621847`)
NEG4

#Insert the vector into a matrix
test1=cbind(NEG4)
test1

# Replicate the column for the number of samples
test2 = test1[,rep(1,7)]
test2

#Subtract values from original matrix using the matrix based on repetition of negative control
PCRbatch4_clean<-PCR_batch_matrix-test2
PCRbatch4_clean

# replace all negative numbers with 0, since we can't have negative seqs reads 
PCRbatch4_cleaned <- replace(PCRbatch4_clean, PCRbatch4_clean < 0, 0)
PCRbatch4_cleaned[,1] # check that it worked

#Now positive controls need to be included 

PCR_batchPOS = subset_samples(juniperus, PCRbatch == "PCR5") 
PCR_batchPOS

# move PCR_batch to a matrix
PCR_batch_matrix = (as(otu_table(PCR_batchPOS), "matrix"))
PCR_batch_matrix

# Now, we can merge all our new OTU tables back together, and reload as a phyloseq object along with the previous taxa table and metedata, since we haven't removed any samples rows or taxa columns.  Once in phyloseq, we can remove columns/OTUs that are now 0.
# check that all objects are in the required format that phyloseq uses for each of them. Otu tables need to be all matrixes, metadata needs to be a data frame and taxonomy table needs to be a matrix 
# now merge otu tables: 
PCR_allbatches_cleaned <- cbind(PCRbatch1_cleaned, PCRbatch2_cleaned, PCRbatch3_cleaned, PCRbatch4_cleaned, PCR_batch_matrix) 

dim(PCR_allbatches_cleaned) # should rows and columns of match what went in

# Create a new phyloseq object with all samples again.

# Export tax table (as matrix) and meta table (as data frame) to be able to re-import them with the new otu table

metadata_matrix = (as(sample_data(juniperus), "matrix"))
metadata_matrix
metadata_df = as.data.frame(metadata_matrix)
rownames(metadata_df)
colnames(metadata_df)

taxdata_matrix = (as(tax_table(juniperus), "matrix"))
taxdata_matrix


# Reimport the clean otu table, meta data and tax table into a new phylseq object
juniperus_cleaned <- phyloseq(otu_table(PCR_allbatches_cleaned, taxa_are_rows=TRUE), 
                             sample_data(metadata_df), 
                             tax_table(taxdata_matrix))

# Check the object to make sure all the samples and OTUs came back.  The columns numbers will be the same, even if they are empty now.
juniperus_cleaned
juniperus

# Clean out OTUs rows that are no longer present.
juniperus_cleaned <- prune_taxa(taxa_sums(juniperus_cleaned) > 0, juniperus_cleaned)

# Clean out samples that are now empty (negative samples)
juniperus_cleaned <- prune_samples(sample_sums(juniperus_cleaned) > 0, juniperus_cleaned)

# Check the objects again
juniperus_cleaned
juniperus

#Check if any OTUs are not present (abundance) in any samples
any(taxa_sums(juniperus_cleaned) == 0)

###################################################################################################### 

# Step 2. DATA TRANSFORMATION

# Transform read counts into  % relative abundances: multiply by 1000 and transform to next integer so it looks like read count
juniperus.rel = transform_sample_counts(juniperus_cleaned, function(x) 100000 * x/sum(x))
otu_table(juniperus.rel) = ceiling(otu_table(juniperus.rel, "matrix")) # transform to next integer so it looks like read count
otu_table(juniperus.rel) # check otu table
juniperus.rel # check project

#check if any OTUs are not present (abundance) in any samples 

any(taxa_sums(juniperus.rel) == 0)

ntaxa(juniperus.rel)

#check distribution of how many reads/OTU, reads/sample 

# Step 3. check distribution of how many reads/OTU, reads/sample: Plot number of reads per OTU / samples 


readsumsdf = data.frame(no.reads = sort(taxa_sums(juniperus), TRUE), sorted = 1:ntaxa(juniperus), type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(no.reads = sort(sample_sums(juniperus), TRUE), sorted = 1:nsamples(juniperus), type = "Samples"))
title = ""
p = ggplot(readsumsdf, aes(x = sorted, y = no.reads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

readsumsdf = data.frame(no.reads = sort(taxa_sums(juniperus.rel), TRUE), sorted = 1:ntaxa(juniperus.rel), type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(no.reads = sort(sample_sums(juniperus.rel), TRUE), sorted = 1:nsamples(juniperus.rel), type = "Samples"))
title = ""
p = ggplot(readsumsdf, aes(x = sorted, y = no.reads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

#Check for distribution of OTUs
hist(log10(taxa_sums(juniperus.rel)))


#################FOR BETA DIVERSITY:  Binomial table: how many OTUS per samples. Filter OTUs from only one sample 

binary_table = transform_sample_counts(juniperus.rel, function(x, minthreshold=0){
  x[x > minthreshold] <- 1
  return(x)})

head(otu_table(binary_table))

#check if any OTUs are not present in any samples 

any(taxa_sums(binary_table) == 0)

ntaxa(binary_table)

#remove OTUs that appear only in 1 sample

any(taxa_sums(binary_table) == 1)

otu_table(prune_taxa(taxa_sums(binary_table) <= 1, binary_table))

subset.1 <- prune_taxa(taxa_sums(binary_table) > 1, binary_table)

subset.1

