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
#SBATCH --time=01:00:00
## Commands to run:

# Load necessary modules
# module load python

# Source .bashrc to ensure conda is initialized
source ~/.bashrc

# Activate the conda environment
conda activate cmd_biopython

# # Debugging: Print conda info and environment list
# conda info
# conda env list

# # Debugging: Print the Python path and check for Biopython
# which python
# python -c "import sys; print(sys.path)"
# python -c "from Bio import SeqIO; print('Biopython is installed and working.')"

# Run the Python script
python contaminate.py

