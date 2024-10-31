#!/bin/bash

# Job name:
#SBATCH --job-name=trim
# Account:
#SBATCH --account=ac_wolflab
# Partition:
#SBATCH --partition=savio2
# Wall clock limit:
#SBATCH --time=00:30:00
# Mail:
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=empchase@berkeley.edu


# Load necessary modules
# module load bio/trimmomatic/0.39-gcc-11.4.0
module load fastqc/0.11.9
module load multiqc/1.17

# # Define the base directory for input, output, and config file
# BASE_DIR="/global/scratch/projects/fc_wolflab/bbioinfo/empchase/synthetic_contaminant"
# # BATCH_DIR="/global/scratch/projects/fc_wolflab/bbioinfo/sarrah/batches"
# # CONFIG_FILE="$BATCH_DIR/slurm_config.csv"

# # head -n 5 $CONFIG_FILE
# # # Extract the line corresponding to this task
# # LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $CONFIG_FILE)
# # echo "$LINE" | awk -F'\t' '{print NF}' 

# # # Parse parameters from the line using a tab as delimiter
# # # Using awk to split and extract fields based on tabs
# # READ1_PATH=$(echo "$LINE" | awk -F'\t' '{print $1}')
# # READ2_PATH=$(echo "$LINE" | awk -F'\t' '{print $2}')
# READ1_PATH="/global/scratch/projects/fc_wolflab/bbioinfo/empchase/synthetic_contaminant/sampleD4+F1_read1.fastq.gz"
# READ2_PATH="/global/scratch/projects/fc_wolflab/bbioinfo/empchase/synthetic_contaminant/sampleD4+F1_read2.fastq.gz"
# # SAMPLE_NAME=$(echo "$LINE" | awk -F'\t' '{print $3}')
# SAMPLE_NAME="sampleD4+F1_90_10"
# # PLATE=$(echo "$LINE" | awk -F'\t' '{print $4}')
# # LIBRARY_NAME=$(echo "$LINE" | awk -F'\t' '{print $5}')
# # ADAPTER_FILE="${BATCH_DIR}/adapter_files/${LIBRARY_NAME}_adapters.fasta"

# # echo "Debug: READ1_PATH: $READ1_PATH"
# # echo "Debug: READ2_PATH: $READ2_PATH"
# # echo "Debug: SAMPLE_NAME: $SAMPLE_NAME"
# # echo "Debug: PLATE: $PLATE"
# # echo "Debug: LIBRARY_NAME: $LIBRARY_NAME"


# # Define output directory
# # OUTPUT_DIR="${BATCH_DIR}/trimmed/${PLATE}/${LIBRARY_NAME}"
# OUTPUT_DIR="trimmomatic_test"


# FASTQC_OUTPUT_DIR="trimmomatic_test"
# # echo "Creating directory at ${OUTPUT_DIR}"

# # mkdir -p "${OUTPUT_DIR}"
# # echo "Creating directory at ${FASTQC_OUTPUT_DIR}"

# # mkdir -p "${FASTQC_OUTPUT_DIR}"


# # Run trimmomatic for paired-end read trimming
# echo "Running Trimmomatic with output in ${OUTPUT_DIR}"
# # echo "Java command length: $(echo -n java -jar /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/trimmomatic/0.36/trimmomatic-0.36.jar PE -phred33 \
# # "$READ1_PATH" \
# # "$READ2_PATH" \
# # "${OUTPUT_DIR}/${SAMPLE_NAME}_output_forward_paired.fq.gz" "${OUTPUT_DIR}/${SAMPLE_NAME}_output_forward_unpaired.fq.gz" \
# # "${OUTPUT_DIR}/${SAMPLE_NAME}_output_reverse_paired.fq.gz" "${OUTPUT_DIR}/${SAMPLE_NAME}_output_reverse_unpaired.fq.gz" \
# # ILLUMINACLIP:"$ADAPTER_FILE":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75 | wc -c)"

# java -jar /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/trimmomatic/0.36/trimmomatic-0.36.jar PE -phred33 \
# "$READ1_PATH" \
# "$READ2_PATH" \
# "${OUTPUT_DIR}/${SAMPLE_NAME}_output_forward_paired.fq.gz" "${OUTPUT_DIR}/${SAMPLE_NAME}_output_forward_unpaired.fq.gz" \
# "${OUTPUT_DIR}/${SAMPLE_NAME}_output_reverse_paired.fq.gz" "${OUTPUT_DIR}/${SAMPLE_NAME}_output_reverse_unpaired.fq.gz" 

# # ILLUMINACLIP:"$ADAPTER_FILE":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75

# # echo "Trimming completed for sample ${SAMPLE_NAME}, Library ${LIBRARY_NAME}"

# #Run fastqc on samples 
# fastqc -o "trimmomatic_test" "${OUTPUT_DIR}/${SAMPLE_NAME}_output_forward_paired.fq.gz"
# fastqc -o "${FASTQC_OUTPUT_DIR}" "${OUTPUT_DIR}/${SAMPLE_NAME}_output_reverse_paired.fq.gz"
# # multiqc -v -o "${FASTQC_OUTPUT_DIR}/${SAMPLE_NAME}_MULTIQC" "${FASTQC_OUTPUT_DIR}"

# echo "Fast QC Completed for sample ${SAMPLE_NAME}, Library ${LIBRARY_NAME}"







# java -jar /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/trimmomatic/0.36/trimmomatic-0.36.jar PE -phred33 \
# /global/scratch/projects/fc_wolflab/bbioinfo/empchase/synthetic_contaminant/sampleD4+F1_read1.fastq.gz \
# /global/scratch/projects/fc_wolflab/bbioinfo/empchase/synthetic_contaminant/sampleD4+F1_read2.fastq.gz \
# trimmomatic_test/sampleD4+f1test_output_forward_paired.fq.gz trimmomatic_test/sampleD4+f1test_output_forward_unpaired.fq.gz \
# trimmomatic_test/sampleD4+f1test_output_reverse_paired.fq.gz trimmomatic_test/sampleD4+f1test_output_reverse_unpaired.fq.gz \
# ILLUMINACLIP:D4+F1_adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75

# #Run fastqc on samples 
fastqc -o trimmomatic_test/fastqc_v0-11-9 trimmomatic_test/sampleD4+f1test_output_forward_paired.fq.gz
fastqc -o trimmomatic_test/fastqc_v0-11-9 trimmomatic_test/sampleD4+f1test_output_reverse_paired.fq.gz
multiqc -v -o trimmomatic_test/fastqc_v0-11-9/multiqc trimmomatic_test/fastqc_v0-11-9