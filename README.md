## Pipeline for analyzing Illumina MiSeq paired-end data from fungal communities

This repository contains the scripts for:

1) Denoising Illumina MiSeq pair-end data (AMPtk)
2) Creating the OTU table and assigning taxonomy (AMPtk)
3) Assigning the trophic mode of OTus (FunGuild inside AMPtk)
4) Filtering the OTU table using R
4) Analyzing alpha- and beta-diversity using R

All the sequence data associated with this project are deposited in [OSF](https://osf.io) while the final OTU table (.biom) is in [data](https://github.com/camillethuyentruong/amptk_pipeline/tree/master/data)


## Citation
Palmer JM, Jusino MA, Banik MT, Lindner DL. 2018. Non-biological synthetic spike-in controls
        and the AMPtk software pipeline improve mycobiome data. PeerJ 6:e4925;
        DOI 10.7717/peerj.4925. https://peerj.com/articles/4925/
