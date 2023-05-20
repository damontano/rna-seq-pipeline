#!/bin/bash

# Define the pipeline name
pipeline_name="RNAgenie Pipeline"

# Display the welcome message
echo "Welcome to $pipeline_name!"

# Set the current date and time
current_date=$(date +"%Y-%m-%d %T")

# Display the current date and time
echo "$pipeline_name started on $current_date"

# Display the system information
echo "System Information:"
uname -a

# Initialize time
SECONDS=0

# Read settings from config file
source config.sh

# Set environment variable for number of threads
export THREADS="$threads"

# Define variables
PRE_TRIMMING_DIR="$working_dir/quality_control/quality_check/pre_trimming"
TRIMMING_DIR="$working_dir/quality_control/trimming"
POST_TRIMMING_DIR="$working_dir/quality_control/quality_check/post_trimming"
MULTIQC_DIR="$working_dir/quality_control/quality_check/multiqc_results"
SAM_DIR="$working_dir/alignment/sam_files"
BAM_DIR="$working_dir/alignment/bam_files"
POST_ALIGNMENT_DIR="$working_dir/quality_control/quality_check/post_alignment"
QUANTIFICATION_DIR="$working_dir/counts"
TMP_DIR="$working_dir/temp"

# Create required directories
mkdir -p "$TMP_DIR" "$PRE_TRIMMING_DIR" "$TRIMMING_DIR" "$POST_TRIMMING_DIR" "$MULTIQC_DIR" "$SAM_DIR" "$BAM_DIR" "$POST_ALIGNMENT_DIR" "$QUANTIFICATION_DIR"

# Check if the user provided input directory
if [ -z "$input_dir" ]; then
	echo "Please provide an input directory with fastq files. Supported file format is .fastq.gz"
	exit 1
fi

# Check if the input directory exists and contains only valid files
if [ ! -d "$input_dir" ]; then
	echo "The input directory $input_dir does not exist."
	exit 1
fi

input_files=$(find "$input_dir" -type f \( -name "*.fastq" \))

if [ -z "$input_files" ]; then
	echo "The input directory $input_dir does not contain any valid fastq files."
	exit 1
fi

# Run RNA-seq pipeline analysis
for FILE in "$input_dir"/*.fastq.gz; do
    SAMPLE_NAME=$(basename "$FILE" .fastq.gz)

    # Check quality of raw reads using FastQC
    fastqc -t $threads -o "$PRE_TRIMMING_DIR" "$FILE"
    echo "FastQC completed for raw reads of sample: $SAMPLE_NAME"

    # Trim sequences using cutadapt
    cutadapt --cores $threads -a $standard_illumina_adapter -m $minimum_read_length -q $quality_cutoff -o "$TRIMMING_DIR"/"$SAMPLE_NAME"_trimmed.fastq.gz "$FILE"
    echo "Trimming completed for sample: $SAMPLE_NAME"

    # Run FastQC on the trimmed sequences
    fastqc -t $threads -o "$POST_TRIMMING_DIR" "$TRIMMING_DIR"/"$SAMPLE_NAME"_trimmed.fastq.gz
    echo "FastQC completed for trimmed sequences of sample: $SAMPLE_NAME"
    
    # Run HISAT2 to align reads
    hisat2 -p $threads -q -x $index_dir -U "$TRIMMING_DIR"/"$SAMPLE_NAME"_trimmed.fastq.gz -S "$SAM_DIR"/"$SAMPLE_NAME"_aligned.sam
    echo "Alignment completed for sample: $SAMPLE_NAME"

    # Convert SAM to sorted BAM for downstream analysis
    samtools view -Sb "$SAM_DIR"/"$SAMPLE_NAME"_aligned.sam | samtools sort -O bam -o "$BAM_DIR"/"$SAMPLE_NAME"_aligned.bam
    #samtools view -Sb "$SAM_DIR"/"$SAMPLE_NAME"_aligned.sam | -o "$BAM_DIR"/"$SAMPLE_NAME"_aligned.bam

    # run FastQC on the aligned sequences
    fastqc -t $threads -o "$POST_ALIGNMENT_DIR" "$SAM_DIR"/"$SAMPLE_NAME"_aligned.sam
    echo "FastQC completed for aligned sequences of sample: $SAMPLE_NAME"

    # Run featureCounts to obtain expression quantification
    featureCounts --tmpDir "$TMP_DIR" -T $threads -a $annotation_dir -o "$QUANTIFICATION_DIR"/count_table.txt "$BAM_DIR"/"$SAMPLE_NAME"_aligned.bam
done

# run MultiQC on all FastQC results
for DIR in "$PRE_TRIMMING_DIR" "$POST_TRIMMING_DIR" "$POST_ALIGNMENT_DIR"; do
    multiqc "$DIR" -o "$MULTIQC_DIR/$(basename "$DIR")"
done

# Display the pipeline completion message and running time
echo "$pipeline_name complete!"
echo "Pipeline finished running in: $((SECONDS / 3600)) hours $((SECONDS / 60 % 60)) minutes $((SECONDS % 60)) seconds."
