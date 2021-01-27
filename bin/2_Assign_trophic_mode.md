# Assigning taxonomy and trophic mode manually using UNITE database

After AMPtk pipeline assignined taxonomy (using [UNITE database](https://unite.ut.ee/)) and fungal guilds (using *FUNGuild*), if there are some OTUs that could not be assigned up to the desired taxonomic level and/or if there were OTUs of interest that were not assigned to a trophic guild or changed in taxonomic resolution, it is possible to perform manual searches to obtain greater resolution. 

The steps to follow to manually add greater taxonomic resolution for your OTUs are indicated in the following lines. 


## Steps 

### Prerequisites

1. Install biom [biom.format.org] (http://biom-format.org/) 

2. If needed, you can add information to the metadata of OTUs and taxonomy in the [biom table](https://github.com/bc-anaisabel/juniperus_paper/blob/master/data/1_taxonomy.biom) obtained in [Script 1](https://github.com/bc-anaisabel/juniperus_paper/blob/master/bin/1_amptk_for_illumina.sh) by using a .txt file like this [one](https://github.com/bc-anaisabel/juniperus_paper/blob/master/data/2_taxonomy.txt). That you can use for the following assignations: 

### Assigning taxonomy 

1. Using the FunGuild table generated with the AMPtk pipeline for ITS2 fungal sequences sort OTUs by most abundant sequence reads and the most frequent (presence/absence across all samples). 

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

2. As the next step is converting to .biom format, and this format requires all taxonomy information in one line, you need to use *concatenate* function in Excel to place all taxonomic ranks and cathegories for each OTU within one cell:

To go from this:
![](../archive/excel_concatenate_2.png)

to this:

![](../archive/excel_concatenate.png)


3. Convert back the file, from *.txt* to *.biom* using:
`biom convert -i taxonomy.txt -o new_taxonomy.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy`

4. If you need more details into how to install the converter consult [Script 3](https://github.com/bc-anaisabel/juniperus_paper/blob/master/bin/3_Convert_txt_to_biom.md) for a full version of instructions. 

5. After this you can work with it in R as follows in [Script 4](https://github.com/bc-anaisabel/juniperus_paper/blob/master/bin/4_Filter_otu_table.R) 

