#!/bin/bash

# Job name:
#SBATCH --job-name=pangenome
# Account:
#SBATCH --account=ac_wolflab
# Partition:
#SBATCH --partition=savio2
# Wall clock limit:
#SBATCH --time=1:00:00
# Mail:
#SBATCH --mail-type=ALL
#SBATCH --mail-user=empchase@berkeley.edu

source ~/.bashrc
conda activate anvio-8


anvi-pan-genome -g PROCHLORO-GENOMES.db \
                --project-name "Prochlorococcus_Pan" \
                --output-dir PROCHLORO \
                --num-threads 6 \
                --minbit 0.5 \
                --mcl-inflation 10 \
                --use-ncbi-blast