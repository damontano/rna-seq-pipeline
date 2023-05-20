# Pipeline for RNA-seq experiments: From Raw Reads to Differential Gene Expression

## Abstract

RNA sequencing (RNA-seq) pipeline is a term used to describe the different steps to analyze gene expression data generated through RNA-seq experiments. The analysis begins with quality check of the raw reads followed by quality control where the raw reads are trimmed. Then, the processing steps include alignment, quantification, normalization, and differential expression analysis. Here we present a simple RNA-seq data analysis pipeline that integrates popular software tools, including FastQC, MultiQC, cutadapt, HISAT2, featurecounts, and DESeq2, to process and analyze RNA-seq data generated from RNA-seq experiments. The pipeline is implemented in a bash script and R Shiny, a web application framework that allows for interactive visualization. We demonstrate the efficacy of our pipeline in analyzing RNA-seq data by identifying differentially expressed genes generated from a mouse mammary gland RNA data. Our pipeline serves as a standardized and efficient platform for high-quality gene expression analysis, which can be further adapted for other RNA-seq datasets.

***

## Project Objectives

1. Create a simple RNA-Seq pipeline to identify differentially expressed genes.
2. Create a public interface to visualize the results of the pipeline.

***

## Background

RNA-seq is a type of high-throughput sequencing technique used to study the transcriptome, which is the complete set of RNA molecules produced by a cell or population of cells. This allows the high-throughput profiling of both coding and non-coding RNA with single-nucleotide resolution. Studying the transcriptome offers researchers the advantage of knowing what genes may or may not be expressed.

While RNA-Seq generates a vast amount of data that can be difficult to analyze, it offers several valuable insights into transcriptomic analyses. Below are some of the most common applications of RNA-Seq data:

1. Changes in gene expression
2. De Novo Transcript Assembly
3. Novel Transcript Discovery

RNA-seq technique is now more accessible to researchers thanks to innovations in next generation sequencing technology, which have significantly decreased the cost of sequencing. Analyzing the RNA-seq data is the final step in the experiment, and choosing the appropriate software is crucial for performing quality control measures, correcting biases, and ultimately answering the research questions posed by the study. The complete steps involved in the RNA workflow are presented in the following figure:

<img src="rna_workflow.png"  alt="Alt text for image"  width="229.96900000000005px"  height="71.305px"  style="object-fit:cover"/>

### Motivation

The motivation for developing the RNA-seq pipeline and R Shiny dashboard is primarily driven by the need to analyze and gain insights into gene expression data. Analyzing RNA sequencing data can be challenging and complicated. Thereby, another motivation for creating the RNA-seq pipeline and R Shiny dashboard to visualize the results includes the ability to generate comprehensive and customizable analyses that are tailored to specific research questions.

### Significance

The significance of creating an RNA-seq pipeline and R Shiny dashboard lies in the ability to generate user-friendly visualizations of RNA-seq data. R Shiny dashboards offer interactive and easy-to-use data visualization to present the results of the RNA-seq pipeline. Therefore, R Shiny dashboards would allow researchers to have visualizations that can present gene expression data information in a more usable, and hopefully  comprehensive way to wider audiences. In addition, the bash pipeline will automate the analysis of RNA sequencing data by incorporating a variety of tools. The goal of the pipeline is to make the analysis of RNA-seq data more efficient and accurate.

***

## Project components

For the project, I decided to use the following components:

1. A reproducible bash shell script.
2. A public facing interface.

***

## Documentation

### Running RNA-seq pipeline

To run the script with the RNA-seq pipeline, the user needs to follow these steps:

1. Download the two .sh extension scripts found in the pipeline folder. The pipeline.sh script contains the entire RNA-seq pipeline. The config.sh script contains the configuration of the pipeline with the appropriate parameters and options for analysis.
2. Open a terminal or command prompt and navigate to the directory that contains the pipeline.sh and the config.sh scripts.
3. Make the scripts executable using the chmod command:

```bash
chmod +x pipeline.sh
chmod +x config.sh
```

4. Run the script using the ./ prefix and the name of the script file:

```bash
./pipeline.sh
```

When the pipeline.sh script runs, it will read the parameters and options from the config.sh script and process the input files as specified. If there are any errors or warnings, the script will display them in the terminal output. When the script finishes, it will display the total time elapsed for the analysis. The current option for mapping and expression quantification is the mouse genome. Users can modify the config.sh script to change input files (must be in fastq.gz format, the index files for alignment, expression quantification.

### Command-line requirements

The pipeline.sh script only accepts data in the format of `fastq.gz`. It also only accepts single-end data that is unstranded.

If the user does not provide an input directory, the following message will be displayed:

```bash
Please provide an input directory with fastq files. Supported file format is .fastq.gz
```

If the input directory does not contain `.fastq.gz` files but the user specified the directory as input file, the following error will be displayed:

```bash
The input directory /path/provided/with/fastq_files does not contain any valid fastq files.
```

If the user sets an input directory that does not exist, the following will be displayed:

```bash
The input directory /path/provided/with/fastq_files does not exist.
```

### Downloading Data

To download RNA-seq data, the SRA toolkit can be used. An example using a for loop on how to download multiple files is shown below:

```bash
for run_id in $(cat /path/to/run_id_list.txt); do fasterq-dump --split-files $run_id --gzip -O fastq_files; done
```
To be able to execute the above command the user must create or have a txt file where the run ids are on a separate line.

### Reference Genome for Alignment

One of the steps of the RNA-seq workflow requires samples to be aligned to a reference genome. For this pipeline, HISAT2 is used to map the reads to the reference genome. HISAT2 provides a list of reference genomes that have already been indexed and can be found [here](http://daehwankimlab.github.io/hisat2/download/). To download the index files provided in HISAT2, the following steps can be performed:

```bash
# Download index for mouse and save it in a new directory
wget -P /path/to/new/directory/ "https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz"
gunzip mm10_genome.tar.gz
tar -xzvf mm10_genome.tar
```
Take note of the name that you give when saving the index to the directory as the path to the directory found in the config.sh script can be updated with the information from the directory. It is also important to note that when specifying the path of the config.sh script for the index file using HISAT2, the name of the base name of the index files must be specified. For example, if we download the mouse index files, and our username is user1, the path in the config.sh script must be specified as follow:

```bash
index_dir=/home/user1/rna_pipeline/genome_index/mm10/genome
```

### Reference Genome for Annotation

When annotating your own data, you will have to provide a reference genome in gtf format. Genomes can be retrieved from the Ensembl FTP site. If you click [here]("https://useast.ensembl.org/index.html"), you will be able to navigate to the Esembl website and find your organism of interest. By right-clicking on the name of the gtf you will be able to copy the URL and then paste into the command-line and download the file. For example, if we were to download the mouse genome we can follow the following steps:

```bash
# Download genome and save it in a new directory
wget -P /path/to/new/directory/ "https://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/Mus_musculus.GRCm39.109.gtf.gz"
gunzip Mus_musculus.GRCm39.109.gtf.gz
```

### Important notes

#### featureCounts

By default, the pipeline.sh script will not count mulit-overlapping reads. This decision was based on the [SubreadUsersGuide.pdf]("https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf"), where for RNA experiments they recommended not to use multi-overlapping.

***

<!---->

## Users

The proposed RNA-seq pipeline using Bash and R-Shiny to visualize the results has the potential to be used by researchers and bioinformaticians. For example, researchers investigating gene expression changes in may find this pipeline useful for analyzing their RNA-seq data. Additionally, bioinformaticians who want to learn how to create an RNA-seq pipeline using Bash and R-shiny may find this pipeline valuable. Overall, the possible users of the pipeline would likely include people who are interested in learning or want to analyze RNA-seq data.
