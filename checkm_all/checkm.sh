#!/bin/bash

# Job name:
#SBATCH --job-name=checkm
# Account:
#SBATCH --account=ac_wolflab
# Partition:
#SBATCH --partition=savio2_bigmem
# Wall clock limit:
#SBATCH --time=24:00:00
# Mail:
#SBATCH --mail-type=ALL
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
input_fasta=/global/scratch/projects/fc_wolflab/bbioinfo/sarrah/batches/trimmed_all/checkm_genomes
output_dir=/global/scratch/projects/fc_wolflab/bbioinfo/empchase/checkm_all/results

# Check if the output directory exists, if so, delete its contents
if [ -d "$output_dir" ]; then
    rm -rf "$output_dir"/*
else
    mkdir -p "$output_dir"
fi

echo $PATH

# Run the CheckM lineage workflow
checkm lineage_wf -x fasta "$input_fasta" "$output_dir" --threads 18 --pplacer_threads 10
