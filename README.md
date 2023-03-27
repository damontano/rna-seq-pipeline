# RNA Seq Pipeline

## Project Overview
This pipeline is designed to take RNA sequencing data from organisms with a reference genome. The RNA seq analysis consists of 5 steps:

Step 1: Visualization of raw reads using FastQC
Step 2: Trimming low quality reads using Trimmomatic
Step 3: Align raw reads to reference genome
Step 4: Count number of reads per annotated gene
Step 5: Analysis of differential gene expression

## Dataset
To test the pipeline, raw data from the NCBI project under the BioProject identifier PRJNA168994 was used. 

## Data visualization
The pipeline runs in the command line where data cannot be visualized. To visualize output files, you can use the provided Shiny app to analyze your data. 

## Requirements
To run this pipeline, the first step is to obtain the raw reads. The following dependencies also need to be installed:
FastQC 
Trimmomatic 
STAR
featureCounts
