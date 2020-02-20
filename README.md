#README

Contains script to process NGS amplicon data using USEARCH and VSEARCH, it can also be used to process any NGS amplicon data. It can handle Ion Torrent, MiSeq, and 454 data. Produces a taxonomy table and a biom table that can be used for downstream analysis with Phyloseq in R. 

First it's required to do a manual installation of amptk using Conda and [USEARCH v9.2.64](http://www.drive5.com/usearch/download.html) from developer.


## Citation
`Palmer JM, Jusino MA, Banik MT, Lindner DL. 2018. Non-biological synthetic spike-in controls
        and the AMPtk software pipeline improve mycobiome data. PeerJ 6:e4925;
        DOI 10.7717/peerj.4925. https://peerj.com/articles/4925/`