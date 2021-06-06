# Data

Here are all the metadata, output and input data files. 

These data files were created and used by the numbered order. 

[1_taxonomy.biom](https://github.com/bc-anaisabel/juniperus_paper/tree/master/data/1_taxonomy.biom) : otu and taxonomy tables obtained from AMPtk pipeline before adding new taxonomy and fungal trophic mode information

[2_taxonomy.txt](https://github.com/bc-anaisabel/juniperus_paper/tree/master/data/2_taxonomy.txt) : otu and taxonomy table in .txt format obtained from AMPtk pipeline that can be edited to add categories of taxonomy and trophic guild for all OTUs 

3_new_tax_biom : this is a residual file, this is not to be used and should be disregarded. I will delete this file once I update the numbers in my scripts. 

[4_new_tax.biom](https://github.com/bc-anaisabel/juniperus_paper/tree/master/data/4_new_tax.biom) : otu and taxonomy tables that have been edited to contain new taxonomy and fungal trophic information and that are now converted back to **.biom**. This is the file to call initially in the first R script: [4_Filter_otu_table.R](https://github.com/bc-anaisabel/juniperus_paper/tree/master/bin/4_Filter_otu_table.R)
   
   

# Metadata 

[amptk.mapping_file.csv](https://github.com/bc-anaisabel/juniperus_paper/tree/master/data/amptk.mapping_file.csv): contains all sample data to be used as input for initial R script. 

The first 13 columns of this file are all different identifiers of each sample. The following 12 columns contain all the information obtained from **soil measurements**: 
- pH; 
- Available phosphorus (mg/kg): Pdis; 
- Calcium (interchangeable cation cmol/kg): Ca;
- Magnesium (interchangeable cation cmol/kg): Mg;
- Potassium (interchangeable cation cmol/kg): K;
- Sodium (interchangeable cation cmol/kg): Na; 
- Hydrogen (interchangeable cation cmol/kg): H; 
- Aluminium (interchangeable cation cmol/kg): Al; 
- Soil moisture (%): SoilM; 
- Carbon (%): C; 
- Nitrogen (%): Nit; 
- Carbon to Nitrogen ratio: CN


[ampt.mapping_file.txt](https://github.com/bc-anaisabel/juniperus_paper/tree/master/data/amptk.mapping_file.txt): contains all sample data to be edited as input for initial R script. This is the same as the amptk.mappin_file.csv file but in .txt format for editing 

[metadata_soil_roots.xlsx](https://github.com/bc-anaisabel/juniperus_paper/tree/master/data/metadata_soil_roots.xlsx): contains soil characterization for soil samples and mycorrhizal colonization information for all samples 
  
  



