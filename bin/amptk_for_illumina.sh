#!/bin/bash

# Script to process demultiplexed Illumina MiSeq 2x300 read, using the amptk pipeline

# amptk needs to be previously installed using conda (https://amptk.readthedocs.io/en/latest/index.html)
### Load the amptk environment with `conda activate amptk`
### Install ITS database using `amptk install -i ITS`



# Pre-processing FASTQ files (forwards and reverse)

amptk illumina -i ../data/fastq -o amptk/ -f GTGARTCATCRARTYTTTG -r CCTSCSCTTANTDATATGC -l 300 --min_len 150 --full_length --cleanup

# --min_len: minimum lenght of sequences = 150 bp ; -l: max lenght = 300
# -f primer forward = gITS7ngs
# -r reverse ITS4ngsUni
# A metadat .txt file is generated at this point and can be completed according to the project


# Move output files to new folder "amptk"

mkdir ../data/amptk
mv amptk.* ../data/amptk


# Clustering at 97% using UPARSE algorithm. This will generate a log file, an OTU table in .txt format and a reference sequence file in fasta

amptk cluster -i ../data/amptk/amptk.demux.fq.gz -o cluster --uchime_ref ITS -m 10


# Move output files to new folder "cluster"

mkdir ../data/cluster
mv cluster* ../data/cluster

# Filter and store in folder
amptk filter -i ../data/cluster/cluster.otu_table.txt -f ../data/cluster/cluster.cluster.otus.fa -o filter -p 0.005 --min_reads_otu 10


# Move output files to new folder "filter"

mkdir ../data/filter
mv filter* ../data/filter


# Assign taxonomy to OTUs

amptk taxonomy -i ../data/filter/filter.final.txt -f ../data/filter/filter.filtered.otus.fa -o taxonomy -m ../data/amptk/amptk.mapping_file.txt -d ITS2 --tax_filter Fungi

# Move output files to new folder "taxonomy"

mkdir ../data/taxonomy
mv taxonomy* ../data/taxonomy

# Assign functional information to OTUs with FunGuild

amptk funguild -i ../data/taxonomy/taxonomy.otu_table.taxonomy.txt -o funguild

# Move output files to new folder "funguild"

mkdir ../data/funguild
mv funguild* ../data/funguild
