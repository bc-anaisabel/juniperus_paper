# Fungal communities in *Juniperus* spp and *Quercus* spp

## Background and sample information ### 

I am working with Illumina MiSeq paired-end data using the ITS2 rDNA region for metabarcoding of fungal communities. In this project I am using metabarcoding to identify the fungal communities present in the soil and in the seedling roots of isolated and mixed oak (*Quercus rugosa*) and juniper (*Juniperus deppeana*) populations. The goal is to study plant-soil feedbacks along a disturbance gradient in central Mexico. 

Because oak and juniper are plants with different mycorrhizal types, oak being known as an ectomycorrhizal plant and juniperus as an arbuscular-mycorrhizal plant, I expect the presence of one influences the other through these plant-soil feedbacks causing changes, which could be positive, neutral, or negative, to the plant growth and overall health and of course, to the fungal communities abundance and diversity. 

The general and bigger outlook is in the context of reforestation where these two plants grow together, I want to know if the fungal communities are affected through plant-soil feedbacks and if this causes any changes that might be relevant to consider when improving reforestation practices and management strategies. Forest management might be improved through studying the biology of these interactions, therefore I consider that in the context of reforestation the study of fungal communities could be a key element for success. 

The samples I used for sequencing were root and soil samples. Samples were taken from 3 different sites: i) disturbed site with a population of *J. deppeana*, ii) mixed site where *J. deppeana* and *Q. rugosa* grow side by side (regeneration zone), and iii) a native site in forest dominated by *Q. rugosa*. Root samples come from collecting the root system of 6 seedlings in each site of *J. deppeana* (sites i and ii) and *Q. rugosa* (sites ii and iii), for a total of 24 root samples. For soil samples, we collected 3 soil cores for each plant species in each site, for a total of 12 soil samples. 

## Repository guide ### 

In this repository you can find scripts, data, metadata and results to identify fungal communities using Illumina MiSeq paired-end data and to analyse fungal diversity and fungal community composition. 

### `/bin`

The [bin](https://github.com/bc-anaisabel/juniperus_paper/tree/master/bin) directory contains the bash, manual steps and R scripts to denoise Illumina MiSeq pair-end data, create the OTU table, assign taxonomy and trophic mode of the OTUs, filter the OTU table in R, and analyze alpha- and beta-diversity. 

#### Scripts and manual steps   

Here you can find the scripts for pre-processing Illumina MiSeq data with AMPtk, the manual steps for assigning OTU taxonomy and trophic mode, and converting metadata from *.txt* to *.biom* to import into R, as well as the numbered (in order of use) R scripts for filtering the OTU table and performing alpha- and beta- diversity analyses. 

The steps to follow are listed in order by letters and within step: **C) R scripts** you can find the R scripts named by numbers in order of use.  

**A) Preprocessing, assigning taxonomy and trophic mode:**

*Preprocessing and taxonomy and trophic mode assignment within AMPtk*

[amptk_for_illumina.sh](https://github.com/bc-anaisabel/juniperus_paper/tree/master/bin/amptk_for_illumina.sh): This is a **bash** script to **denoise** Illumina MiSeq pair-end data, **create an OTU table**, and **assign taxonomy** and fungal **trophic guilds** within **AMPtk**. [AMPtk_pipeline.md](https://github.com/bc-anaisabel/juniperus_paper/tree/master/bin/AMPtk_pipeline.md): this is a text file that describes what the *amptk_for_illumina.sh* script does. 


**B) Manually completing taxonomy, trophic mode and converting your file after running AMPtk pipeline**

[Assign_trophic_mode.md](https://github.com/bc-anaisabel/juniperus_paper/tree/master/bin/Assign_trophic_mode.md): instructions for manually assigning taxonomy and trophic mode to OTUs that could not be identified beyond Order (or higher up in taxonomy rank) in AMPtk using UNITE database.
  
If you need to edit your file to add the new taxonomic and trophic assigments, you will need to go from a .txt file to a .biom file:
  
[Convert_txt_to_biom.md](https://github.com/bc-anaisabel/juniperus_paper/tree/master/bin/Convert_txt_to_biom.md): instructions for manually converting the edited **[2_taxonomy.txt](https://github.com/bc-anaisabel/juniperus_paper/tree/master/data/2_taxonomy.txt)** into **[4_new_tax_biom](https://github.com/bc-anaisabel/juniperus_paper/tree/master/data/4_new_tax.biom)**, the file that you will be able to use in R with *phyloseq* package. 

**C) R scripts** 

Use in number order to run alpha- and beta- diversity analyses as follows:
  * [1_Filter_otu_table.R](https://github.com/bc-anaisabel/juniperus_paper/tree/master/bin/1_Filter_otu_table.R) : filtering otu table
  * [2_Juniperus_Alpha_Diversity.R](https://github.com/bc-anaisabel/juniperus_paper/tree/master/bin/2_Juniperus_Alpha_Diversity.R): performing alpha diversity analysis
  * [3_Juniperus_Beta_Diversity.R](https://github.com/bc-anaisabel/juniperus_paper/tree/master/bin/3_Juniperus_Beta_Diversity.R): performing beta diversity analysis
  * [4_Models.R](https://github.com/bc-anaisabel/juniperus_paper/tree/master/bin/4_Models.R): modelling diversity data 
  

### `/data`

The [data](https://github.com/bc-anaisabel/juniperus_paper/tree/master/data) directory contains the metadata file containing samples information and the output data from the AMPtk pipeline that acts as input data for R scripts. These data files were created and used by the numbered order. 

#### Data files


[1_taxonomy.biom](https://github.com/bc-anaisabel/juniperus_paper/tree/master/data/1_taxonomy.biom) : otu and taxonomy tables obtained from AMPtk pipeline before adding new taxonomy and fungal trophic mode information

[2_taxonomy.txt](https://github.com/bc-anaisabel/juniperus_paper/tree/master/data/2_taxonomy.tx) : otu and taxonomy table in .txt format obtained from AMPtk pipeline that can be edited to add categories of taxonomy and trophic guild for all OTUs 

3_new_tax_biom : this is a residual file, this is not to be used and should be disregarded. I will delete this file once I update the numbers in my scripts. 

[4_new_tax_biom](https://github.com/bc-anaisabel/juniperus_paper/tree/master/data/4_new_tax_biom) : otu and taxonomy tables that have been edited to contain new taxonomy and fungal trophic information and that are now converted back to **.biom**. This is the file to call initially in the first R script: [1_Filter_otu_table.R](https://github.com/bc-anaisabel/juniperus_paper/tree/master/bin/1_Filter_otu_table.R)
   
   

#### Metadata 

amptk.mapping_file.csv: contains all sample data to be used as input for initial R script 

ampt.mapping_file.txt: contains all sample data to be edited as input for initial R script 

metadata_soil_roots.xlsx: contains soil characterization for soil samples and mycorrhizal colonization information for all samples 


### `/output`

The [output](https://github.com/bc-anaisabel/juniperus_paper/tree/master/output) directory contains all figures obtained from R scripts to be used for the publication. 


### Raw data ### 


All the sequence data associated with this project are deposited in [OSF](https://osf.io) while the final OTU table (.biom) is in [data](https://github.com/bc-anaisabel/juniperus_paper/tree/master/data)


### Citation
Palmer JM, Jusino MA, Banik MT, Lindner DL. 2018. Non-biological synthetic spike-in controls
        and the AMPtk software pipeline improve mycobiome data. PeerJ 6:e4925;
        DOI 10.7717/peerj.4925. https://peerj.com/articles/4925/
