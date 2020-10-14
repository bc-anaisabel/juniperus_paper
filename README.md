## Fungal communities in *Juniperus spp* and *Quercus spp* 

I am working with Illumina MiSeq paired-end data using the ITS2 rDNA region for metabarcoding of fungal communities. 

In this repository you can find scripts for:

1. Identify fungal communities ussing Illumina MiSeq paired-end data 

* 1.1 Denoising Illumina MiSeq pair-end data (AMPtk)
* 1.2 Creating the OTU table and assigning taxonomy (AMPtk)
* 1.3 Assigning the trophic mode of OTUs (FunGuild inside AMPtk)

2. Analyse fungal diversity and fungal community composition 

* 2.1 Filtering the OTU table using R
* 2.2 Analyzing alpha- and beta-diversity using R

All the sequence data associated with this project are deposited in [OSF](https://osf.io) while the final OTU table (.biom) is in [data](https://github.com/bc-anaisabel/juniperus_paper/tree/master/data)


## Citation
Palmer JM, Jusino MA, Banik MT, Lindner DL. 2018. Non-biological synthetic spike-in controls
        and the AMPtk software pipeline improve mycobiome data. PeerJ 6:e4925;
        DOI 10.7717/peerj.4925. https://peerj.com/articles/4925/
