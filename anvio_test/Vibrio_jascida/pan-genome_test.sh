#!/bin/bash

# Job name:
#SBATCH --job-name=EC_checkmtest
# Account:
#SBATCH --account=ac_wolflab
# Partition:
#SBATCH --partition=savio2
# Wall clock limit:
#SBATCH --time=00:20:00
# Mail:
#SBATCH --mail-type=ALL
#SBATCH --mail-user=empchase@berkeley.edu

#BEFORE RUNNING: run in terminal: conda activate anvio-8-2

# Code to run

anvi-pan-genome -g V_jascida_genomes/V_jascida-GENOMES.db \
                --project-name V_jascida \
                --num-threads 4