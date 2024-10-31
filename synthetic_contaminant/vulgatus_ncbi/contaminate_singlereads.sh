#!/usr/bin/env bash
#Author: empchase@berkeley.edu

#SBATCH --job-name=contaminate
#
# Account:
#SBATCH --account=ac_wolflab
# Partition:
#SBATCH --partition=savio2
#
# Wall Clock Limit:
#SBATCH --time=02:00:00
## Commands to run:

# # Activate the conda environment
# conda activate cmd_biopython

# Define file names for python script
fasta1_sampleA_read1=/global/scratch/projects/fc_wolflab/bbioinfo/empchase/synthetic_contaminant/vulgatus_ncbi/vulgatus_ATCC8482_pseudoreads.fasta
fasta1_sampleB_read1=/global/scratch/projects/fc_wolflab/bbioinfo/empchase/synthetic_contaminant/vulgatus_ncbi/vulgatus_CL09T03C04_pseudoreads.fasta

output_subsetA_read1="subsetATCC8482_pseudo.fasta"
output_subsetB_read1="subsetCL09T03C04_pseudo.fasta"

proportion_A=0.9 # also define fraction of first file that will go into contamination 

# final output filenames for bash concatentation
output_combined_read1="ATCC8482+CL09T03C04_pseudo.fasta"




# # Run the Python script
python contaminate_singlereads.py \
    $fasta1_sampleA_read1 \
    $output_subsetA_read1 \
    --proportion_A $proportion_A \
    $fasta1_sampleB_read1 \
    $output_subsetB_read1 \

# Concatenate read1 files from Sample A and Sample B
cat "$output_subsetA_read1" "$output_subsetB_read1" > "$output_combined_read1"

echo "Concatenation complete. Output files: $output_combined_read1 and $output_combined_read2"

