
Here are all the metadata and output and input data files. 

=================
## Data


**1_taxonomy.biom** : otu and taxonomy tables obtained from AMPtk pipeline before adding new taxonomy and fungal trophic mode information

**2_taxonomy.txt** : otu and taxonomy table in .txt format obtained from AMPtk pipeline that has been edited to add categories of 
   taxonomy and trophic mode of all OTUs 

**3_new_tax_biom** : delete 

**4_new_tax_biom** : otu and taxonomy tables that have been edited to contain new taxonomy and fungal trophic information and converted back to .biom
   This is the file to call initially in * 1_Filter_otu_table.R *
   
   

## Metadata 

**amptk.mapping_file.csv**: contains all sample data to be used as input for initial R script ( same as *amptk.mapping_file.txt* ) 

**metadata_soil_roots.xlsx**: contains soil characterization for soil samples and mycorrhizal colonization information for all samples 



## Text files in order of use 

**1_Assign_trophic_mode.md**: text file of instructions for assigning taxonomy and trophic mode to OTUs that could not be identified beyond Order 
  (or higher up in taxonomy rank) using UNITE database. 

**2_Convert_txt_to_biom.md**: text file of instructions for converting *2_taxonomy.txt* into **.biom** format to use in R with the *phyloseq* package
