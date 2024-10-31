#!/usr/bin/env bash
#Author: empchase@berkeley.edu

#SBATCH --job-name=pseduoreads
#
# Account:
#SBATCH --account=ac_wolflab
# Partition:
#SBATCH --partition=savio2
#
# Wall Clock Limit:
#SBATCH --time=00:10:00
## Commands to run:

# Variables for input, output, read length, and step size
INPUT_GENOME="GCF_000273295.1_Bact_vulg_CL09T03C04_V1_genomic.fna.gz"
OUTPUT_PSEUDO_READS="vulgatus_CL09T03C04_pseudoreads.fasta"
READ_LENGTH=100
STEP_SIZE=75

# Call the Python script
python pseudoreads.py "$INPUT_GENOME" "$OUTPUT_PSEUDO_READS" --read_length $READ_LENGTH --step_size $STEP_SIZE
