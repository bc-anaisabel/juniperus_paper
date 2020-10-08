# Assigning taxonomy and trophic mode manually using UNITE database



1. Using the FunGuild table generated with the AMPtk pipeline for ITS2 fungal sequences (**Step 6**) sort OTUs by most abundant sequence reads and the most frequent (presence/absence across all samples). 

2. For the top 100 OTUs that were not identified beyond Order (including those that were not identified taxonomically at Phylum or Class level), search FASTA sequences in UNITE database doing a blast search for each OTU sequence individually: 
 
* Parameters: 

**blast e-value**: 1 (default)

**Datasets to include**: INSD and environmental 

**Check**: ID %, SH (**97**%) 

3. Provide SH and assign taxonomy when possible with criteria:

* 90% similarity: Family, 95% Genus, 97% Species

## Assigning trophic mode 

For all OTUs, when FunGuild did not assign a trophic mode or if taxonomy was updated to a finer level with a new SH: 

1. Check SH (Species Hypothesis) at 97% and assign trophic mode using UNITE Search Page (Ecology tab). If there is no SH associated use the INSD accession number. If UNITE does not provide a trophic mode or status as either mycorrhizal or not mycorrhizal, leave as **unknown**.  

## Update metadata (e.g. taxonomy, trophic mode, soil details) in .biom table

### Install biom [biom.format.org] (http://biom-format.org/) 

1. Convert to a .txt file( *taxonomy.txt* )

2. Edit in Excel with all previous steps. Create one column with the mycorrhizal/no mycorrhizal cathegory and one with the trophic mode. 

3. Biom format requires all taxonomy information as one line. Use concatenate function in Excel to place all taxonomic ranks and cathegories within one cell. 

3. Convert back to biom format to work with it in R  

