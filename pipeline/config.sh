#!/bin/bash

# paths
working_dir=rna_seq_pipeline
input_dir="{PATH_TO_FASTQ_FILES}"
index_dir="{PATH_TO_INDEX_GENOME}"
annotation_dir="{PATH_TO_ANNOTATION_GTF}"

# adapters
standard_illumina_adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC

# cutadapt_parameters
minimum_read_length=20
quality_cutoff=20

# parameters
threads=8
