# Scripts and manual steps  

Here you can find the scripts for pre-processing Illumina MiSeq data with AMPtk, the manual steps for assigning OTU taxonomy and trophic mode, and converting metadata from *.txt* to *.biom* to import into R, as well as the numbered (in order of use) R scripts for filtering the OTU table and performing alpha- and beta- diversity analyses. 

The steps to follow are listed in order:  

- ## Step 1: Preprocessing, assigning taxonomy and trophic mode:


[1_amptk_for_illumina.sh](https://github.com/bc-anaisabel/juniperus_paper/tree/master/bin/1_amptk_for_illumina.sh): This is a **bash** script to **denoise** Illumina MiSeq pair-end data, **create an OTU table**, and **assign taxonomy** and fungal **trophic guilds** within **AMPtk**. [AMPtk_pipeline.md](https://github.com/bc-anaisabel/juniperus_paper/tree/master/bin/AMPtk_pipeline.md): this is a text file that describes what the *1_amptk_for_illumina.sh* script does. 

- ## Step 2: Manually completing taxonomy, trophic mode and converting your file after running AMPtk pipeline 

[2_Assign_trophic_mode.md](https://github.com/bc-anaisabel/juniperus_paper/tree/master/bin/2_Assign_trophic_mode.md): instructions for manually assigning taxonomy and trophic mode to OTUs that could not be identified beyond Order (or higher up in taxonomy rank) by AMPtk with the UNITE database.
  
- ## Step 3: If you edited your file as in Step 2, you will need to change your edited *.txt* file to a *.biom* file:
  
[3_Convert_txt_to_biom.md](https://github.com/bc-anaisabel/juniperus_paper/tree/master/bin/3_Convert_txt_to_biom.md): instructions for manually converting the edited **[2_taxonomy.txt](https://github.com/bc-anaisabel/juniperus_paper/tree/master/data/2_taxonomy.txt)** into **[4_new_tax.biom](https://github.com/bc-anaisabel/juniperus_paper/tree/master/data/4_new_tax.biom)**, the file that you will be able to use in R with *phyloseq* package. 

In my case I edited added a two more categories to the taxonomy line, so in my[2_taxonomy.txt](https://github.com/bc-anaisabel/juniperus_paper/tree/master/data/2_taxonomy.txt) file I used the following abbreviations to denote mycorrhizal category (t) and trophic category(a): 

- t__myc = mycorrhizal fungi;
	
	a__ecm = ectomycorrhizal; 
	
	a__am = arbuscular mycorrhizal;
	
	a__otro = mycorrhizal fungi other than ectomycorrhizal or arbuscular;
	
	a__unknown = mycorrhizal but not identified to a mycorrhizal category


- t__nomyc = non-mycorrhizal fungi;
	
	a__sap = saprophytic;
	
	a__endo = endophytic;
	
	a__par = parasitic; 
	
	a__lic = lichen;
	
	a__unknown = unknown category of non-mycorrhizal fungi
	
- t__unknown = fungi that could not be identified either as mycorrhizal or non-mycorrhizal 

Once I converted the *.txt* file to the biom file: [4_new_tax.biom](https://github.com/bc-anaisabel/juniperus_paper/tree/master/data/4_new_tax.biom), all of my taxonomic information included the listed fungal categories. So throught the following **R scripts** you'll run into those abbreviations when I am trying to group and plot specific groups of fungi. 

- ## Step 4: R scripts 
Use in number order to run alpha- and beta- diversity analyses as follows:
  * 4_Filter_otu_table.R : filtering otu table
  * 5_Juniperus_Alpha_Diversity.R: performing alpha diversity analysis
  * 6_Juniperus_Beta_Diversity.R: performing beta diversity analysis
  * 7_networks.R: modelling diversity data 
  
 
  
  



  

