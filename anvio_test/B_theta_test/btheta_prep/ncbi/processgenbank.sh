#!/bin/bash

# Job name:
#SBATCH --job-name=process_genbank
# Account:
#SBATCH --account=ac_wolflab
# Partition:
#SBATCH --partition=savio2
# Wall clock limit:
#SBATCH --time=00:30:00
# Mail:
#SBATCH --mail-type=ALL
#SBATCH --mail-user=empchase@berkeley.edu
source ~/.bashrc
conda activate anvio-8

anvi-script-process-genbank -i DSM_2079.gbff \
                            --output-gene-calls DSM_2079_gene_calls.tsv \
                            --output-functions DSM_2079_functions.tsv \
                            --output-fasta DSM_2079_refs.fa \
                            --annotation-source prodigal

for genome in $(ls *.fna | cut -f1 -d ".")
do

  anvi-script-reformat-fasta -o "$genome"_clean.fa \
                             --prefix "$genome" \
                             --simplify-names \
                             "$genome".fna
  prodigal -f gff \
           -c \
           -i "$genome"_clean.fa \
           -o "$genome".gff
done