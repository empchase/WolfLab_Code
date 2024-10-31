#!/bin/bash

# Job name:
#SBATCH --job-name=anvi_btheta
# Account:
#SBATCH --account=ac_wolflab
# Partition:
#SBATCH --partition=savio2
# Wall clock limit:
#SBATCH --time=00:10:00
# Mail:
#SBATCH --mail-type=ALL
#SBATCH --mail-user=empchase@berkeley.edu

#make sure anvio-8-2 conda env is activated
anvi-import-misc-data binding_shigella_extrainfo.txt \
                              -p B_theta/B_theta-PAN.db \
                              --target-data-table layers

anvi-compute-functional-enrichment-in-pan -p B_theta/B_theta-PAN.db \
                                          -g B_theta-GENOMES.db \
                                          --category Galactan_Binder \
                                          --annotation-source KOfam \
                                          -o binding_enriched-functions.txt

anvi-compute-functional-enrichment-in-pan -p B_theta/B_theta-PAN.db \
                                          -g B_theta-GENOMES.db \
                                          --category Shigella_response \
                                          --annotation-source COG20_FUNCTION \
                                          -o shigella_enriched-functions.txt                                          