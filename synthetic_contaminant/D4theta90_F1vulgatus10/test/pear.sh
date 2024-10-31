#!/usr/bin/env bash
#Author: empchase@berkeley.edu

#SBATCH --job-name=pear
#
# Account:
#SBATCH --account=ac_wolflab
# Partition:
#SBATCH --partition=savio2
#
# Wall Clock Limit:
#SBATCH --time=02:00:00
## Commands to run:

# # Load necessary modules
# module load pear

# pear -f /global/scratch/projects/fc_wolflab/bbioinfo/empchase/synthetic_contaminant/subsetD4_read1.fastq.gz -r /global/scratch/projects/fc_wolflab/bbioinfo/empchase/synthetic_contaminant/subsetD4_read2.fastq.gz -o test_outputs/subsetD4_2
pear -f /global/scratch/projects/fc_wolflab/bbioinfo/230914_BacteroidesPatnodeLab_GTAC_VPL/BacteroidesPatnodePlate1_121323/LIB015053_22GMJHLT3_S57_L006_R1_001.fastq.gz -r /global/scratch/projects/fc_wolflab/bbioinfo/230914_BacteroidesPatnodeLab_GTAC_VPL/BacteroidesPatnodePlate1_121323/LIB015053_22GMJHLT3_S57_L006_R2_001.fastq.gz -o test_outputs/LIB015053_22GMJHLT3_S57_L006



# readone=$1
# readtwo=$2
# result=$3

# module load pear

# pear -f $readone -r $readtwo -o $result -p 0.001