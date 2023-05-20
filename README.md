# RNA Seq Pipeline

## Project Overview
This pipeline is designed to take RNA sequencing data from organisms with a reference genome. The RNA seq analysis consists of 5 steps:

- Step 1: Quality check of raw reads using FastQC
- Step 2: Trimming low quality reads using Cutadapt
- Step 3: Align raw reads to reference genome
- Step 4: Count number of reads per annotated gene
- Step 5: Analysis of differential gene expression using R-shiny

## Dataset
To test the pipeline, raw data from the Gene Expression Omnibus database (GEO) under accession number GSE60450 was used.

## Data visualization
The pipeline runs in the command line where data cannot be visualized. To visualize output files, you can upload the count table obtained as the last step in the pipeline to perform differential expression analysis. 

## Requirements
To run this pipeline, the first step is to obtain the raw reads in fastq format. Then, save the fastq files in your home directory in a folder name fastq_files. In addition to the fasta files, the following dependencies need to be installed:
- FastQC 
- Trimmomatic 
- featureCounts
- featureCounts
