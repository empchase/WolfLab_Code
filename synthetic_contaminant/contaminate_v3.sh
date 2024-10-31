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
fastq1_sampleA_read1="/global/scratch/projects/fc_wolflab/bbioinfo/230914_BacteroidesPatnodeLab_GTAC_VPL/BacteroidesPatnodePlate1_121323/LIB015015_22GMJHLT3_S95_L006_R1_001.fastq.gz"
fastq2_sampleA_read2="/global/scratch/projects/fc_wolflab/bbioinfo/230914_BacteroidesPatnodeLab_GTAC_VPL/BacteroidesPatnodePlate1_121323/LIB015015_22GMJHLT3_S95_L006_R2_001.fastq.gz"
fastq1_sampleB_read1="/global/scratch/projects/fc_wolflab/bbioinfo/230914_BacteroidesPatnodeLab_GTAC_VPL/BacteroidesPatnodePlate1_121323/LIB015018_22GMJHLT3_S92_L006_R1_001.fastq.gz"
fastq2_sampleB_read2="/global/scratch/projects/fc_wolflab/bbioinfo/230914_BacteroidesPatnodeLab_GTAC_VPL/BacteroidesPatnodePlate1_121323/LIB015018_22GMJHLT3_S92_L006_R2_001.fastq.gz"

output_subsetA_read1="subsetA02_read1.fastq.gz"
output_subsetA_read2="subsetA02_read2.fastq.gz"
output_subsetB_read1="subsetA05_read1.fastq.gz"
output_subsetB_read2="subsetA05_read2.fastq.gz"

proportion_A=0.9 # also define fraction of first file that will go into contamination 

# final output filenames for bash concatentation
output_combined_read1="sampleA02+A05_read1.fastq.gz"
output_combined_read2="sampleA02+A05_read2.fastq.gz"



# # Run the Python script
python contaminate_v3.py \
    $fastq1_sampleA_read1 \
    $fastq2_sampleA_read2 \
    $output_subsetA_read1 \
    $output_subsetA_read2 \
    --proportion_A $proportion_A \
    $fastq1_sampleB_read1 \
    $fastq2_sampleB_read2 \
    $output_subsetB_read1 \
    $output_subsetB_read2

# Concatenate read1 files from Sample A and Sample B
cat "$output_subsetA_read1" "$output_subsetB_read1" > "$output_combined_read1"

# Concatenate read2 files from Sample A and Sample B
cat "$output_subsetA_read2" "$output_subsetB_read2" > "$output_combined_read2"

echo "Concatenation complete. Output files: $output_combined_read1 and $output_combined_read2"

