#!/bin/bash

# Create directory
mkdir fastq_files
cd fastq_files

# Download short version of fastq files as specified in
# https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/rna-seq-reads-to-counts/tutorial.html#mapping
for i in `seq 2444 2455`; do
    wget "https://zenodo.org/record/4249555/files/SRR155$i.fastq.gz"
done
