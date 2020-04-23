# AMPtk pipeline for ITS2 fungal sequences

The script *amptk_for_illumina.sh* runs the following comands:


## Step 1. 

AMPtk needs to be previously installed [using conda](https://amptk.readthedocs.io/en/latest/index.html)

The usearch9 dependency must also be manually installed from the [developer](http://www.drive5.com/usearch/download.html) and linked to the PATH.

Load the amptk environment with: `conda activate amptk`

Install ITS database using: `amptk install -i ITS`
 

## Step 2. 
### Pre-processing FASTQ files (forward and reverse): 

amptk illumina <arguments>: 

-i, --fastq Input folder of FASTQ files (Required)

-o, --out Ouput folder name. Default: amptk-data

-f, --fwd_primer Forward primer sequence. Default: fITS7

-r, --rev_primer Reverse primer sequence. Default: ITS4

--min_len Minimum length read to keep. Default: 100

--full_length Keep only full length sequences

--cleanup Remove intermediate files

`amptk illumina -i ../data/fastq -o amptk/ -f GTGARTCATCRARTYTTTG -r CCTSCSCTTANTDATATGC -l 300 --min_len 150 --full_length --cleanup`

*Minimum lenght of sequences = 150 bp ; max lenght = 300, primer forward = gITS7ngs, reverse ITS4ngsUni*


### Move output files to new folder "amptk":

`mkdir ../data/amptk`

`mv amptk.* ../data/amptk`


## Step 3. 
### Clustering at 97% using UPARSE algorithm: 

* The first step is to quality trim the reads using expected-errors quality trimming. The quality trimmed reads are then pushed through each of the algorithms which output a set of OTUs. The method used for clustering here is UPARSE. This will generate a log file, an OTU table in .txt format and a reference sequence file in fasta

amptk cluster <arguments>: 

-i, --fastq Input FASTQ file (Required)

-o, --out Ouput base name. Default: out

-m, --minsize Minimum number of reads per OTU (singleton filter). Default: 2

--uchime_ref Run Ref Chimera filtering. Default: off [ITS, LSU, COI, 16S, custom path]

`amptk cluster -i ../data/amptk/amptk.demux.fq.gz -o cluster --uchime_ref ITS -m 10`


### Move output files to new folder "cluster":

`mkdir ../data/cluster`

`mv cluster* ../data/cluster`

## Step 4. 
### Filter OTU table to combat barcode-switching or index-bleed:

* Index-bleed (reads assigned to the wrong sample) seems to happen during the sequencing process of Illumina. AMPtk uses spike-in-artificial mock communities to measure the degree of index-bleed in a sequencing run and then conservatively applies that threshold to remove read counts that are within the range of index-bleed from an OTU table. If a spike-in mock was not included in the run, this threshold can be defined manually or by default (0.005 % of reads assigned to the wrong sample)

amptk filter <arguments>:

-i, --otu_table OTU table generated from the amptk cluster

-f, --fasta OTU fasta generated from the amptk cluster

-o, --out Base name for output files. Defaults: use input basename

-p, --index_bleed Filter index bleed between samples (percent). Default: 0.005

--min_reads_otu Minimum number of reads for valid OTU to be retained. Default:2


`amptk filter -i ../data/cluster/cluster.otu_table.txt -f ../data/cluster/cluster.cluster.otus.fa -o filter -p 0.005 --min_reads_otu 10`


### Move output files to new folder "filter": 

`mkdir ../data/filter`

`mv filter* ../data/filter`



## Step 5.
### Assign taxonomy to OTUs:

* The default method in AMPtk is a 'hybrid' taxonomy assignment which calculates a consensus LCA (last common ancestor) taxonomy based on the results of Global Alignment (USEARCH/VSEARCH), UTAX, and SINTAX. The LCA method used here is applied twice, first in finding the LCA of identical Global Alignment hits and then again to compare the taxonomy from Global Alignment to the Bayesian Classifiers. If there is no conflict between taxonomies from each method, the taxonomy string with most information is retained (SINTAX/UTAX results are used if BLAST-like search pct identity is less than 97%%. If % identity is greater than 97%, the result with most taxonomy levels is retained).

AMPtk uses the [UNITE database](https://unite.ut.ee/) for ITS taxonomy assignement. Since UNITE is a curated database, the taxonomy assignment usullay gives much better results than GenBank (used for example by QIMME).

amptk taxonomy <arguments>:

-i, --otu_table Input OTU table file (i.e. otu_table from amptk cluster)

-f, --fasta Input FASTA file (i.e. OTUs from amptk cluster). (Required)

-o, --out Base name for output file. Default: amptk-taxonomy.<method>.txt
 
-m, --mapping_file QIIme-like mapping file

-d, --db Select pre-installed database [ITS1, ITS2, ITS, 16S LSU, COI]. Default: ITS2

--tax_filter Remove OTUs from OTU table that do not match filter, i.e. Fungi to keep only fungi


`amptk taxonomy -i ../data/filter/filter.final.txt -f ../data/filter/filter.filtered.otus.fa -o taxonomy -m ../data/amptk/amptk.mapping_file.txt -d ITS2 --tax_filter Fungi`

*At this step, the mapping file geneated in the `amptk illumina` (Step 2) that has been completed with meta data information can be called to include the information in the final output .biom file

### Move output files to new folder "taxonomy":

`mkdir ../data/taxonomy`

`mv taxonomy* ../data/taxonomy`

## Step 6. 
### Assign functional information to OTUs with FunGuild:

* For this step an OTU table that has been appended with taxonomy needs to be provided.

amptk funguild <arguments>:

-i, --input Input OTU table

-d, --db Database to use [fungi, nematode]. Default: fungi

-o, --out Output file basename

`amptk funguild -i ../data/taxonomy/taxonomy.otu_table.taxonomy.txt -o funguild`

### Move output files to new folder "funguild":

`mkdir ../data/funguild`

`mv funguild* ../data/funguild`

## Step 7. 
### Run analyses with R using .biom files generated from "taxonomy" folder

Filtering of the OTU table (negative controls, OTUs and samples) are realized in R (see script filter_biom.r in bin).
