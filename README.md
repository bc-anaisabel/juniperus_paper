## Fungal communities in *Juniperus* spp and *Quercus* spp

### Background and sample information ### 

I am working with Illumina MiSeq paired-end data using the ITS2 rDNA region for metabarcoding of fungal communities. In this project I am using metabarcoding to identify the fungal communities present in the soil and in the seedling roots of isolated and mixed oak (*Quercus rugosa*) and juniper (*Juniperus deppeana*) populations. The goal is to study plant-soil feedbacks along a disturbance gradient in central Mexico. 

Because oak and juniper are plants with different mycorrhizal types, oak being known as an ectomycorrhizal plant and juniperus as an arbuscular-mycorrhizal plant, it is expected that the presence of both may influence the other through these plant-soil feedbacks causing changes, which could be positive, neutral, or negative, to the plant growth and overall health and of course, to the fungal communities abundance and diversity. 

The samples I used for sequencing were root and soil samples. Samples were taken from 3 different sites: i) disturbed site with a population of *J. deppeana*, ii) mixed site where *J. deppeana* and *Q. rugosa* grow side by side (regeneration zone), and iii) a native site in forest dominated by *Q. rugosa*. Root samples come from collecting the root system of 6 seedlings in each site of *J. deppeana* (sites i and ii) and *Q. rugosa* (sites ii and iii), for a total of 24 root samples. For soil samples, we collected 3 soil cores for each plant species in each site, for a total of 12 soil samples. 

### Repository guide ### 

In this repository you can find scripts, data, metadata and results to identify fungal communities using Illumina MiSeq paired-end data and to analyse fungal diversity and fungal community composition. 

The `/bin` directory contains the bash, manual steps and R scripts to denoise Illumina MiSeq pair-end data, create the OTU table, assign taxonomy and trophic mode of the OTUs, filter the OTU table in R, and analyze alpha- and beta-diversity. 

The `/data` directory contains the metadata file containing samples information and the output data from the AMPtk pipeline that acts as input data for R scripts. It also contains the .md files numbered in order of use to transform the output data into input data (from *.biom* to *.txt* and back to *.biom*).  

The `/output` directory contains all figures obtained from R scripts to be used for the publication. 

### Raw data ### 


All the sequence data associated with this project are deposited in [OSF](https://osf.io) while the final OTU table (.biom) is in [data](https://github.com/bc-anaisabel/juniperus_paper/tree/master/data)


### Citation
Palmer JM, Jusino MA, Banik MT, Lindner DL. 2018. Non-biological synthetic spike-in controls
        and the AMPtk software pipeline improve mycobiome data. PeerJ 6:e4925;
        DOI 10.7717/peerj.4925. https://peerj.com/articles/4925/
