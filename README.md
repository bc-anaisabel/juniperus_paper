## Fungal communities in *Juniperus* spp and *Quercus* spp

I am working with Illumina MiSeq paired-end data using the ITS2 rDNA region for metabarcoding of fungal communities. 

In this repository you can find scripts, data, metadata and results to identify fungal communities using Illumina MiSeq paired-end data and to analyse fungal diversity and fungal community composition. 

The `/bin` directory contains the bash and R scripts for: 

1. Denoising Illumina MiSeq pair-end data (AMPtk)
2. Creating the OTU table and assigning taxonomy (AMPtk)
3. Assigning the trophic mode of OTUs (FunGuild inside AMPtk)
4. Filtering the OTU table using R
5. Analyzing alpha- and beta-diversity using R

The `/data` directory contains the metadata file that gathers all sample information and the output data from the AMPtk pipeline that acts as input data for R scripts. It also contains the .md files numbered in order of use to transform the output data into input data.  

The `/output` directory contains all figures obtained from R scripts to be used for the publication. 



All the sequence data associated with this project are deposited in [OSF](https://osf.io) while the final OTU table (.biom) is in [data](https://github.com/bc-anaisabel/juniperus_paper/tree/master/data)


## Citation
Palmer JM, Jusino MA, Banik MT, Lindner DL. 2018. Non-biological synthetic spike-in controls
        and the AMPtk software pipeline improve mycobiome data. PeerJ 6:e4925;
        DOI 10.7717/peerj.4925. https://peerj.com/articles/4925/
