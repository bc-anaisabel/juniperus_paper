# Scripts 
===========

## Scripts 

Scripts for pre-processing Illumina MiSeq data with AMPtk, assigning OTU taxonomy and trophic mode, updating metadata from .txt to .biom to import into R, as well as the numbered R scripts for filtering OTU table and performing alpha- and beta- diversity analyses. 


1. *amptk_for_illumina.sh* is the first-step. This is a bash script to denoise Illumina MiSeq pair-end data, create an OTU table, and assign taxonomy and fungal trophic modes within AMPtk

2. R scripts 1-4 can be used by number order to run alpha- and beta- diversity analyses as follows:
  * 1_Filter_otu_table.R : 
  * 2_Juniperus_Alpha_Diversity.R:
  * 3_Juniperus_Beta_Diversity.R:
  * 4_Models.R
  
  ## Text files 
  
  *AMPtk_pipeline.md* is the text file that describes what the *amptk_for_illumina.sh* script does. 
  
  

