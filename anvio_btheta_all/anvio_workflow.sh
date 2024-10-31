#!/bin/bash

# Job name:
#SBATCH --job-name=anvi_btheta
# Account:
#SBATCH --account=ac_wolflab
# Partition:
#SBATCH --partition=savio2
# Wall clock limit:
#SBATCH --time=10:00:00
# Mail:
#SBATCH --mail-type=ALL
#SBATCH --mail-user=empchase@berkeley.edu

# download the reference
# someday here: a curl command, but for now I'm download+scp ; download into the same folder as your other fasta files

# Initiate variables
SPECIES=B_theta
# REFERENCE_SCAFFOLDS=GCF_014131755.1_ASM1413175v1_genomic.fna
# # rename reference genome
# mv $REFERENCE_SCAFFOLDS ref_scaffolds.fasta

# # you MUST run the reformat fasta on a reference genome,
# #for loop / prep
# ls *.fasta | awk 'BEGIN{FS="."}{print $1}' > genomes.txt
# for g in `cat genomes.txt`
# do
#     echo
#     echo "Working on $g ..."
#     echo
#     anvi-script-reformat-fasta ${g}.fasta \
#                                --min-len 1000 \
#                                --simplify-names \
#                                -o ${g}_reformat.fasta
# done

# # generate contigs db
# for g in `cat genomes.txt`
# do
#     echo
#     echo "Working on $g ..."
#     echo
#     anvi-gen-contigs-database -f ${g}_reformat.fasta \
#                               -o ${SPECIES}_${g}.db \
#                               --num-threads 18 \
#                               -n ${SPECIES}_${g}
# done

# #Annotate contigs db

# # Comment out if you've already run these?
# anvi-setup-ncbi-cogs --num-threads 18
# anvi-setup-scg-taxonomy --num-threads 18

# # annotate
# for g in *.db
# do
#     anvi-run-hmms -c $g --num-threads 18
#     anvi-run-ncbi-cogs -c $g --num-threads 18 
#     anvi-scan-trnas -c $g --num-threads 18
#     anvi-run-scg-taxonomy -c $g --num-threads 18 
#     anvi-estimate-scg-taxonomy -c $g --debug
# done

# # generate text file describing where your databases are, used by anvio
# anvi-script-gen-genomes-file --input-dir . -o external-genomes.txt
# # check completeness for contamination
# anvi-estimate-genome-completeness -e external-genomes.txt

# generate genome storage db (the prep db for pangenome analysis)
anvi-gen-genomes-storage -e external-genomes.txt \
                         -o ${SPECIES}-GENOMES.db


anvi-pan-genome -g ${SPECIES}-GENOMES.db \
    --project-name ${SPECIES} \
    --num-threads 18 

anvi-compute-genome-similarity --external-genomes external-genomes.txt \
                               --program pyANI --output-dir ANI \
                               --num-threads 18 \
                               --pan-db ${SPECIES}/${SPECIES}-PAN.db


