#!/bin/bash

# Job name:
#SBATCH --job-name=EC_prokkatest
# Account:
#SBATCH --account=ac_wolflab
# Partition:
#SBATCH --partition=savio2
# Wall clock limit:
#SBATCH --time=1:00:00
# Mail:
#SBATCH --mail-type=ALL
#SBATCH --mail-user=empchase@berkeley.edu


module load python
source activate biopython

echo "biopython activated"

prokka --prefix EC_test /global/scratch/projects/fc_wolflab/bbioinfo/sarrah/batches/trimmed/Plate1/A1/75/assembly/contigs.fasta