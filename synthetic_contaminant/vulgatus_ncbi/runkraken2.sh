#!/bin/bash
# Job name:
#SBATCH --job-name=krake
# Account:
#SBATCH --account=ac_wolflab
# Partition:
#SBATCH --partition=savio3_bigmem

#SBATCH --time=01:00:00
# Mail:
#SBATCH --mail-type=END
#SBATCH --mail-user=empchase@berkeley.edu


#Define stuff
KRAKENDB="/global/scratch/projects/fc_wolflab/software/kraken2/kraken_prebuilt_db"

READ1_PATH=/global/scratch/projects/fc_wolflab/bbioinfo/empchase/synthetic_contaminant/vulgatus_ncbi/ATCC8482+CL09T03C04_pseudo.fasta
SAMPLE_NAME=ATCC8482+CL09T03C04_pseudo

echo "Processing sample ${SAMPLE_NAME} with reads ${READ1_PATH} "

# Run Kraken2 and log the command
kraken2 --db ${KRAKENDB} --threads 16 --minimum-hit-groups 3 \
--report "$SAMPLE_NAME".k2report \
$READ1_PATH  > ${SAMPLE_NAME}.kraken2


