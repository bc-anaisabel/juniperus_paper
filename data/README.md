# Data

Here are all the metadata, output and input data files. 

These data files were created and used by the numbered order. 


**1_taxonomy.biom** : otu and taxonomy tables obtained from AMPtk pipeline before adding new taxonomy and fungal trophic mode information

**2_taxonomy.txt** : otu and taxonomy table in .txt format obtained from AMPtk pipeline that has been edited to add categories of 
   taxonomy and trophic mode of all OTUs 

**3_new_tax_biom** : this is a residual file, this is not to be used and should be disregarded. I will delete this file once I update the numbers in my scripts. 

**4_new_tax_biom** : otu and taxonomy tables that have been edited to contain new taxonomy and fungal trophic information and converted back to .biom
   This is the file to call initially in the R script located in `/bin` [1_Filter_otu_table.R](https://github.com/bc-anaisabel/juniperus_paper/tree/master/bin/1_Filter_otu_table.R)
   
   

## Metadata 

**amptk.mapping_file.csv**: contains all sample data to be used as input for initial R script 

**ampt.mapping_file.txt**: contains all sample data to be edited as input for initial R script 

**metadata_soil_roots.xlsx**: contains soil characterization for soil samples and mycorrhizal colonization information for all samples 



