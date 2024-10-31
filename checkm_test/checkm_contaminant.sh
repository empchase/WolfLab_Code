#!/bin/bash

# Job name:
#SBATCH --job-name=checkm_contam
# Account:
#SBATCH --account=ac_wolflab
# Partition:
#SBATCH --partition=savio2
# Wall clock limit:
#SBATCH --time=00:20:00
# Mail:
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=empchase@berkeley.edu

export PATH=$PATH:/global/home/users/empchase/.local/bin
export MODULEPATH=${MODULEPATH}:/clusterfs/vector/home/groups/software/sl-7.x86_64/modfiles
export CHECKM_DATA_PATH=/global/scratch/projects/fc_wolflab/software/checkm_data

# # Load the CheckM and HMMER modules
# module load checkm
# module load hmmer
module load prodigal
module load pplacer/1.1a19

# Define the input directory and output directory
input_fastadir="/global/scratch/projects/fc_wolflab/bbioinfo/empchase/synthetic_contaminant/A8thetaVPI-90_A4thetaWH502-10"
output_dir="/global/scratch/projects/fc_wolflab/bbioinfo/empchase/checkm_test/contamBtheta_A8vA4_results_rocky8"

# Check if the output directory exists, if so, delete its contents
if [ -d "$output_dir" ]; then
    rm -rf "$output_dir"/*
else
    mkdir -p "$output_dir"
fi

# Run the CheckM lineage workflow
checkm lineage_wf -x fasta "$input_fastadir" "$output_dir"

# # Run the CheckM quality assessment
# checkm qa "$output_dir/lineage.ms" "$output_dir" #future can get rid of this line bc checkm lineage_wf runs qa

