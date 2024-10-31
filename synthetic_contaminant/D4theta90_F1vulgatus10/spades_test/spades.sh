#!/usr/bin/env bash
#Author: empchase@berkeley.edu

#SBATCH --job-name=spades
#
# Account:
#SBATCH --account=ac_wolflab
# Partition:
#SBATCH --partition=savio2
#
# Wall Clock Limit:
#SBATCH --time=01:00:00
# Mail:
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=empchase@berkeley.edu
# Commands to run:
# Source .bashrc to ensure conda is initialized
source ~/.bashrc
# Activate the conda environment
conda activate python3-8-8
export MODULEPATH=$MODULEPATH:/global/scratch/projects/fc_wolflab/software/modfiles/

module load spades/3.13.0
spades.py  -m 46 -o . --pe1-1 /global/scratch/projects/fc_wolflab/bbioinfo/empchase/synthetic_contaminant/trimmomatic_test/trimmomatic_test/sampleD4+f1test_output_forward_paired.fq.gz --pe1-2 /global/scratch/projects/fc_wolflab/bbioinfo/empchase/synthetic_contaminant/trimmomatic_test/trimmomatic_test/sampleD4+f1test_output_reverse_paired.fq.gz
