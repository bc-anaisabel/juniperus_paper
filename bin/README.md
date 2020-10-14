# Scripts 
===========

Scripts for pre-processing Illumina MiSeq data with AMPtk, assigning OTU taxonomy and trophic mode, updating metadata from .txt to .biom to import into R, as well as the numbered R scripts for filtering OTU table and performing alpha- and beta- diversity analyses. 

##Scripts 

1. *amptk _for _illumina.sh* is the initial script to denoising Illumina MiSeq pair-end data, create an OTU table, and assign taxonomy and fungal trophic modes within AMPtk

2. OTUs that could not be identified beyond Order (or higher up in taxonomy rank) can be individually blasted using UNITE database following instructions from **Assign _trophic _mode.md**. Additionally, missing trophic modes can also be checked with these instructions in UNITE. 

3. Metadata in .txt format can be edited and turned back into .biom to use in phyloseq in R following **Convert_ biom_ to_ txt.md**

4. R scripts 1-4 can be used by number order to run alpha- and beta- diversity analyses 

