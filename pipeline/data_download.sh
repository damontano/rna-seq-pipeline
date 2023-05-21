#!/bin/bash

# Create directory
mkdir fastq_files
cd fastq_files

# Download short version of fastq files as specified in
for i in `seq 2444 2455`; do
    wget "https://zenodo.org/record/4249555/files/SRR155$i.fastq.gz"
done