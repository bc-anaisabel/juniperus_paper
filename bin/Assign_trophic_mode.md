# Assigning taxonomy and trophic mode manually using UNITE database

After AMPtk pipeline assignined taxonomy (using *UNITE* database) and fungal guilds (using *FUNGuild*), if there are some OTUs that could not be assigned up to the desired taxonomic level and/or if there were OTUs of interest that were not assigned to a trophic guild or changed in taxonomic resolution, it is possible to perform manual searches to obtain greater resolution. The steps to follow in order to obtain this greater taxonomic resolution are indicated in the following lines. 


## Steps 

### Prerequisites

1. Install biom [biom.format.org] (http://biom-format.org/) 

2. If needed, update metadata in .biom table obtained in xxxxxxxx 

3. Convert .biom (insert here name of biom file) to a .txt file ( *taxonomy.txt* )

### Assigning taxonomy 

1. Using the FunGuild table generated with the AMPtk pipeline for ITS2 fungal sequences (**Step 6**) sort OTUs by most abundant sequence reads and the most frequent (presence/absence across all samples). 

2. For the top 100 OTUs that were not identified beyond Order (including those that were not identified taxonomically at Phylum or Class level), search FASTA sequences in UNITE database doing a blast search for each OTU sequence individually: 
 
* Parameters: 

**blast e-value**: 1 (default)

**Datasets to include**: INSD and environmental 

**Check**: ID %, SH (**97**%) 

3. Provide SH and assign taxonomy when possible with criteria:

* 90% similarity: Family, 95% Genus, 97% Species


### Assigning trophic mode 

For all OTUs, when FunGuild did not assign a trophic mode or if taxonomy was updated to a finer level with a new SH: 

1. Check SH (Species Hypothesis) at 97% and assign trophic mode using UNITE Search Page (Ecology tab). If there is no SH associated use the INSD accession number. If UNITE does not provide a trophic mode or status as either mycorrhizal or not mycorrhizal, leave as **unknown**. 


### Assembling new taxonomic and fungal guilds information 

1. Update .txt file in Excel with the new taxonomic and fungal guild information obtained. You can choose whatever format you need to do this, for example I added one column with the mycorrhizal/no mycorrhizal cathegory and one with the trophic mode. 

2. As the next step is converting to .biom format, and this format requires all taxonomy information in one line, you need to use *concatenate* function in Excel to place all taxonomic ranks and cathegories for each OTU within one cell. 

3. Convert back the file, from .txt to .biom and now you can work with it in R  

