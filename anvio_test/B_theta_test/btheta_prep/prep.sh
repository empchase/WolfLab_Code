#!/bin/bash

# Job name:
#SBATCH --job-name=pangenome
# Account:
#SBATCH --account=ac_wolflab
# Partition:
#SBATCH --partition=savio2
# Wall clock limit:
#SBATCH --time=00:20:00
# Mail:
#SBATCH --mail-type=ALL
#SBATCH --mail-user=empchase@berkeley.edu

source ~/.bashrc
conda activate anvio-8

contig_path=/global/scratch/projects/fc_wolflab/bbioinfo/sarrah/batches/trimmed_all/Plate1/C4/B_theta_WH507_assembly/contigs.fasta
sample=C4_WH507




anvi-script-reformat-fasta $contig_path -o ${sample}_assembly_contigs-fasta --simplify-names --prefix $sample   

singularity exec /global/scratch/projects/fc_wolflab/bbioinfo/prokka.sif prokka --prefix PROKKA_${sample} --outdir PROKKA_${sample} --cpus 2 ${sample}_assembly_contigs-fasta 

python /global/scratch/projects/fc_wolflab/software/gff_parser.py PROKKA_${sample}/PROKKA_${sample}.gff --gene-calls ${sample}_gene_calls.txt --annotation ${sample}_gene_annot.txt  

anvi-gen-contigs-database -f ${sample}_assembly_contigs-fasta -o ../${sample}_contigs.db --external-gene-calls ${sample}_gene_calls.txt --project-name ${sample}_ContigDB --ignore-internal-stop-codons


