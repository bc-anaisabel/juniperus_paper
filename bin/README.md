# Scripts and manual steps  

Here you can find the scripts for pre-processing Illumina MiSeq data with AMPtk, the manual steps for assigning OTU taxonomy and trophic mode, and converting metadata from *.txt* to *.biom* to import into R, as well as the numbered (in order of use) R scripts for filtering the OTU table and performing alpha- and beta- diversity analyses. 

The steps to follow are listed in order by letters and within step: **C) R scripts** you can find the R scripts named by numbers in order of use.  

## A) Preprocessing, assigning taxonomy and trophic mode:

**Preprocessing and taxonomy and trophic mode assignment within AMPtk**

[amptk_for_illumina.sh](https://github.com/bc-anaisabel/juniperus_paper/tree/master/bin/amptk_for_illumina.sh): This is a **bash** script to **denoise** Illumina MiSeq pair-end data, **create an OTU table**, and **assign taxonomy** and fungal **trophic guilds** within **AMPtk**. [AMPtk_pipeline.md](https://github.com/bc-anaisabel/juniperus_paper/tree/master/bin/AMPtk_pipeline.md): this is a text file that describes what the *amptk_for_illumina.sh* script does. 

## B) Manually completing taxonomy, trophic mode and converting your file after running AMPtk pipeline** 

[Assign_trophic_mode.md](https://github.com/bc-anaisabel/juniperus_paper/tree/master/bin/Assign_trophic_mode.md): instructions for manually assigning taxonomy and trophic mode to OTUs that could not be identified beyond Order (or higher up in taxonomy rank) in AMPtk using UNITE database.
  
If you need to edit your file to add the new taxonomic and trophic assigments, you will need to go from a .txt file to a .biom file:
  
[Convert_txt_to_biom.md](https://github.com/bc-anaisabel/juniperus_paper/tree/master/bin/Convert_txt_to_biom.md): instructions for manually converting the edited **[2_taxonomy.txt](https://github.com/bc-anaisabel/juniperus_paper/tree/master/data/2_taxonomy.txt)** into **[4_new_tax_biom](https://github.com/bc-anaisabel/juniperus_paper/tree/master/data/4_new_tax.biom)**, the file that you will be able to use in R with *phyloseq* package. 

## C) R scripts 
Use in number order to run alpha- and beta- diversity analyses as follows:
  * 1_Filter_otu_table.R : filtering otu table
  * 2_Juniperus_Alpha_Diversity.R: performing alpha diversity analysis
  * 3_Juniperus_Beta_Diversity.R: performing beta diversity analysis
  * 4_Models.R: modelling diversity data 
  
 
  
  



  

