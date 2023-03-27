#!/bin/bash

# Create working directory
mkdir rna_seq_pipeline

# Move into directory
cd rna_seq_pipeline

# Step 0: Obtain Data
mv ~/fastq_files ~/rna_seq_pipeline

# Step 1: Quality Control 
# Step 1a: Check quality of reads using FastQC
# Step 1b: Create output folder
mkdir fastqc_output 

# Step 1c: Run FastQC
fastqc fastq_files/*.fastq -o fastqc_output

# Step 1d: Run MultiQC
cd fastq_output
multiqc .


