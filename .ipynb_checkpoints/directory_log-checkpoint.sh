# June 14, 2024
# test out prokka on one of sarrah's assemblies
touch prokka_test.sh
    module load python
    source activate biopython

    echo "biopython activated"

    prokka --prefix EC_test /global/scratch/projects/fc_wolflab/bbioinfo/sarrah/batches/trimmed/Plate1/A1/75/assembly/contigs.fasta
    # running successfully, set the time to one hour, we'll see how it goes


# June 17, 2024
# mess around with prokka outputs

#realized that my directory was called EC, which is the same abbreviation for Enzyme Commision numbers, so I changed to empchase
# from /global/scratch/projects/fc_wolflab/bbioinfo, run:
mv EC empchase

#

# June 27, 2024
# try out checkm

#first organize
mkdir prokka_test
mv prokka_test.sh prokka_test
mv slurm-19342843.out prokka_test/
mv EC_test prokka_test

mkdir checkm_test
touch checkm_test/checkm_test.sh

#can't remember if I did this before but this is where the checkm function lives so :
export MODULEPATH=${MODULEPATH}:home/groups/software/sl-7.x86_64/modfiles

#try checkm
module load checkm #it ran? did it work? not sure
chmod u+rwx checkm_test/checkm_test.sh

# July 1, 2024
# fix checkm script

#last error message:
# OSError: [Errno 2] No such file or directory: '/global/scratch/projects/fc_wolflab/bbioinfo/empchase/checkm_test/checkm_results/bins'
# add a command to add bin if it doesn't exist

# July 2, 2024
# test anvio on a genome

mkdir anvio_test
anvi-gen-contigs-database -f ../prokka_test/EC_test/EC_test.fna -o EC_test_fna.db -n "EC fna test" 

Input FASTA file .............................: /global/scratch/projects/fc_wolflab/bbioinfo/empchase/prokka_test/EC_test/EC_test.fna
Name .........................................: EC fna test
Description ..................................: No description is given
Num threads for gene calling .................: 1                               

Finding ORFs in contigs
===============================================
Genes ........................................: /tmp/tmpvp_i5f8z/contigs.genes
Amino acid sequences .........................: /tmp/tmpvp_i5f8z/contigs.amino_acid_sequences
Log file .....................................: /tmp/tmpvp_i5f8z/00_log.txt

CITATION
===============================================
Anvi'o will use 'prodigal' by Hyatt et al (doi:10.1186/1471-2105-11-119) to
identify open reading frames in your data. When you publish your findings,
please do not forget to properly credit their work.

Result .......................................: Prodigal (v2.6.3) has identified 5992 genes.

                                                                                
CONTIGS DB CREATE REPORT
===============================================
Split Length .................................: 20,000
K-mer size ...................................: 4
Skip gene calling? ...........................: False
External gene calls provided? ................: False
Ignoring internal stop codons? ...............: False
Splitting pays attention to gene calls? ......: True
Contigs with at least one gene call ..........: 585 of 615 (95.1%)              
Contigs database .............................: A new database, EC_test_fna.db, has been created.
Number of contigs ............................: 615
Number of splits .............................: 864
Total number of nucleotides ..................: 7,175,274
Gene calling step skipped ....................: False
Splits broke genes (non-mindful mode) ........: False
Desired split length (what the user wanted) ..: 20,000
Average split length (what anvi'o gave back) .: 21,667

anvi-gen-contigs-database took 0:01:04.229633

# so it gives me a database EC_test_fna.db (expected)'
# run hmms
anvi-run-hmms -c EC_test_fna.db

Contigs DB ...................................: EC_test_fna.db                  
HMM sources ..................................: Ribosomal_RNA_23S, Ribosomal_RNA_16S, Bacteria_71, Archaea_76, Ribosomal_RNA_5S, Ribosomal_RNA_12S, Ribosomal_RNA_18S, Protista_83, Ribosomal_RNA_28S
Alphabet/context target found ................: RNA:CONTIG
Alphabet/context target found ................: AA:GENE                         
Target sequences determined ..................: 615 sequences for RNA:CONTIG; 5,992 sequences for AA:GENE

HMM Profiling for Ribosomal_RNA_23S
===============================================
Reference ....................................: Seeman T, https://github.com/tseemann/barrnap
Kind .........................................: Ribosomal_RNA_23S
Alphabet .....................................: RNA
Context ......................................: CONTIG
Domain .......................................: N/A
HMM model path ...............................: /tmp/tmpw3_j3zs8/Ribosomal_RNA_23S.hmm
Number of genes in HMM model .................: 2
Noise cutoff term(s) .........................: --cut_ga
Number of CPUs will be used for search .......: 1
HMMer program used for search ................: nhmmscan
Temporary work dir ...........................: /tmp/tmp3v94gi55
Log file for thread 0 ........................: /tmp/tmp3v94gi55/RNA_contig_sequences.fa.0_log
                                                                                
Done with Ribosomal_RNA_23S ????

Number of raw hits in table file .............: 2                               
Number of weak hits removed by HMMER parser ..: 0
Number of hits in annotation dict  ...........: 2
Pruned .......................................: 1 out of 2 hits were removed due to redundancy
Gene calls added to db .......................: 1 (from source "Ribosomal_RNA_23S")

HMM Profiling for Ribosomal_RNA_16S
===============================================
Reference ....................................: Seeman T, https://github.com/tseemann/barrnap
Kind .........................................: Ribosomal_RNA_16S
Alphabet .....................................: RNA
Context ......................................: CONTIG
Domain .......................................: N/A
HMM model path ...............................: /tmp/tmpw3_j3zs8/Ribosomal_RNA_16S.hmm
Number of genes in HMM model .................: 3
Noise cutoff term(s) .........................: --cut_ga
Number of CPUs will be used for search .......: 1
HMMer program used for search ................: nhmmscan
Temporary work dir ...........................: /tmp/tmp3v94gi55
Log file for thread 0 ........................: /tmp/tmp3v94gi55/RNA_contig_sequences.fa.0_log
                                                                                
Done with Ribosomal_RNA_16S ????

Number of raw hits in table file .............: 1                               
Number of weak hits removed by HMMER parser ..: 0
Number of hits in annotation dict  ...........: 1
Gene calls added to db .......................: 1 (from source "Ribosomal_RNA_16S")
                                                                                
HMM Profiling for Bacteria_71
===============================================
Reference ....................................: Lee modified, https://doi.org/10.1093/bioinformatics/btz188
Kind .........................................: singlecopy
Alphabet .....................................: AA
Context ......................................: GENE
Domain .......................................: bacteria
HMM model path ...............................: /tmp/tmpw3_j3zs8/Bacteria_71.hmm
Number of genes in HMM model .................: 71
Noise cutoff term(s) .........................: --cut_ga
Number of CPUs will be used for search .......: 1
HMMer program used for search ................: hmmscan
Temporary work dir ...........................: /tmp/tmpb2zaq29c
Log file for thread 0 ........................: /tmp/tmpb2zaq29c/AA_gene_sequences.fa.0_log
                                                                                
Done with Bacteria_71 ????

Number of raw hits in table file .............: 77                              
Number of weak hits removed by HMMER parser ..: 0
Number of hits in annotation dict  ...........: 77
                                                                                
HMM Profiling for Archaea_76
===============================================
Reference ....................................: Lee, https://doi.org/10.1093/bioinformatics/btz188
Kind .........................................: singlecopy
Alphabet .....................................: AA
Context ......................................: GENE
Domain .......................................: archaea
HMM model path ...............................: /tmp/tmpw3_j3zs8/Archaea_76.hmm
Number of genes in HMM model .................: 76
Noise cutoff term(s) .........................: --cut_ga
Number of CPUs will be used for search .......: 1
HMMer program used for search ................: hmmscan
Temporary work dir ...........................: /tmp/tmpb2zaq29c
Log file for thread 0 ........................: /tmp/tmpb2zaq29c/AA_gene_sequences.fa.0_log
                                                                                
Done with Archaea_76 ????

Number of raw hits in table file .............: 40                              
Number of weak hits removed by HMMER parser ..: 0
Number of hits in annotation dict  ...........: 40
                                                                                
HMM Profiling for Ribosomal_RNA_5S
===============================================
Reference ....................................: Seeman T, https://github.com/tseemann/barrnap
Kind .........................................: Ribosomal_RNA_5S
Alphabet .....................................: RNA
Context ......................................: CONTIG
Domain .......................................: N/A
HMM model path ...............................: /tmp/tmpw3_j3zs8/Ribosomal_RNA_5S.hmm
Number of genes in HMM model .................: 5
Noise cutoff term(s) .........................: --cut_ga
Number of CPUs will be used for search .......: 1
HMMer program used for search ................: nhmmscan
Temporary work dir ...........................: /tmp/tmp3v94gi55
Log file for thread 0 ........................: /tmp/tmp3v94gi55/RNA_contig_sequences.fa.0_log
                                                                                
Done with Ribosomal_RNA_5S ????

Number of raw hits in table file .............: 0                               

* The HMM source 'Ribosomal_RNA_5S' returned 0 hits. SAD (but it's OK).'
                                                                                
HMM Profiling for Ribosomal_RNA_12S
===============================================
Reference ....................................: Seeman T, https://github.com/tseemann/barrnap
Kind .........................................: Ribosomal_RNA_12S
Alphabet .....................................: RNA
Context ......................................: CONTIG
Domain .......................................: N/A
HMM model path ...............................: /tmp/tmpw3_j3zs8/Ribosomal_RNA_12S.hmm
Number of genes in HMM model .................: 1
Noise cutoff term(s) .........................: --cut_ga
Number of CPUs will be used for search .......: 1
HMMer program used for search ................: nhmmscan
Temporary work dir ...........................: /tmp/tmp3v94gi55
Log file for thread 0 ........................: /tmp/tmp3v94gi55/RNA_contig_sequences.fa.0_log
                                                                                
Done with Ribosomal_RNA_12S ????

Number of raw hits in table file .............: 0                               

* The HMM source 'Ribosomal_RNA_12S' returned 0 hits. SAD (but it's OK).
                                                                                
HMM Profiling for Ribosomal_RNA_18S
===============================================
Reference ....................................: Seeman T, https://github.com/tseemann/barrnap
Kind .........................................: Ribosomal_RNA_18S
Alphabet .....................................: RNA
Context ......................................: CONTIG
Domain .......................................: N/A
HMM model path ...............................: /tmp/tmpw3_j3zs8/Ribosomal_RNA_18S.hmm
Number of genes in HMM model .................: 1
Noise cutoff term(s) .........................: --cut_ga
Number of CPUs will be used for search .......: 1
HMMer program used for search ................: nhmmscan
Temporary work dir ...........................: /tmp/tmp3v94gi55
Log file for thread 0 ........................: /tmp/tmp3v94gi55/RNA_contig_sequences.fa.0_log
                                                                                
Done with Ribosomal_RNA_18S ????

Number of raw hits in table file .............: 0                               

* The HMM source 'Ribosomal_RNA_18S' returned 0 hits. SAD (but it's OK).
                                                                                
HMM Profiling for Protista_83
===============================================
Reference ....................................: Delmont, http://merenlab.org/delmont-euk-scgs
Kind .........................................: singlecopy
Alphabet .....................................: AA
Context ......................................: GENE
Domain .......................................: eukarya
HMM model path ...............................: /tmp/tmpw3_j3zs8/Protista_83.hmm
Number of genes in HMM model .................: 83
Noise cutoff term(s) .........................: -E 1e-25
Number of CPUs will be used for search .......: 1
HMMer program used for search ................: hmmscan
Temporary work dir ...........................: /tmp/tmpb2zaq29c
Log file for thread 0 ........................: /tmp/tmpb2zaq29c/AA_gene_sequences.fa.0_log
                                                                                
Done with Protista_83 ????

Number of raw hits in table file .............: 2                               
Number of weak hits removed by HMMER parser ..: 0
Number of hits in annotation dict  ...........: 2
                                                                                
HMM Profiling for Ribosomal_RNA_28S
===============================================
Reference ....................................: Seeman T, https://github.com/tseemann/barrnap
Kind .........................................: Ribosomal_RNA_28S
Alphabet .....................................: RNA
Context ......................................: CONTIG
Domain .......................................: N/A
HMM model path ...............................: /tmp/tmpw3_j3zs8/Ribosomal_RNA_28S.hmm
Number of genes in HMM model .................: 1
Noise cutoff term(s) .........................: --cut_ga
Number of CPUs will be used for search .......: 1
HMMer program used for search ................: nhmmscan
Temporary work dir ...........................: /tmp/tmp3v94gi55
Log file for thread 0 ........................: /tmp/tmp3v94gi55/RNA_contig_sequences.fa.0_log
                                                                                
Done with Ribosomal_RNA_28S ????

Number of raw hits in table file .............: 0                               

* The HMM source 'Ribosomal_RNA_28S' returned 0 hits. SAD (but it's OK).
                                                                                
??? anvi-run-hmms took 0:01:37.312058

[M 27


# check the databases
(anvio-8) [empchase@n0027 anvio_test]$ ls
EC_test_fna.db
# so there's still just one db, the last command must have annotated the db

# try pangenome
(anvio-8) [empchase@n0027 anvio_test]$ anvi-gen-genomes-storage -e EC_test_fna.db -o genomes-storage.db


File/Path Error: The file at 'EC_test_fna.db' does not seem to be plain a text file :/


(anvio-8) [empchase@n0027 anvio_test]$ anvi-gen-genomes-storage -e ../prokka_test/EC_test/EC_test.fna -o EC_test_genomes_storage.db -n "EC fna test" 
usage: anvi-gen-genomes-storage [-h] [-e FILE_PATH] [-i FILE_PATH]
                                [--gene-caller GENE-CALLER] -o GENOMES_STORAGE
                                [--version] [--debug] [--force]
                                [--fix-sad-tables] [--quiet] [--no-progress]
                                [--as-markdown] [--display-db-calls]
                                [--force-use-my-tree] [--force-overwrite]
                                [--tmp-dir TMP_DIR]
                                [--I-know-this-is-not-a-good-idea]
anvi-gen-genomes-storage: error: unrecognized arguments: -n EC fna test

(anvio-8) [empchase@n0027 anvio_test]$ anvi-gen-genomes-storage -e ../prokka_test/EC_test/EC_test.fna -o EC_test_genomes_storage.db                  


File/Path Error: You (or some code on your behalf) asked anvi'o if the file at                   
                 '../prokka_test/EC_test/EC_test.fna' was a TAB-delimited file. Anvi'o took the  
                 very first line in it that did not start with the character '#' (as in commented
                 out lines), and found zero TAB in it. This is not how we make TAB-delimited     
                 files :(                                                                        


# checkm retry - going to try giving it a directory as an input instead of a file
*******************************************************************************
 [CheckM - tree] Placing bins in reference genome tree.
*******************************************************************************

  [Error] Output directory must be empty: /global/scratch/projects/fc_wolflab/bbioinfo/empchase/checkm_test/checkm_results


  Controlled exit resulting from an unrecoverable error or warning.

# deleted checkm_results/bins/, trying again

# July 3, 2024
# get checkm to work

# my script was still adding /bins so I deleted that line of code
# latest error (/global/scratch/projects/fc_wolflab/bbioinfo/empchase/checkm_test/slurm-19731388.out) said:
  [Error] Make sure HMMER executables (e.g., hmmsearch, hmmfetch) are on your system path.
#so I added module load hmmer into the script

# July 15, 2024
# create synthetically contaminated file
zcat 230914_BacteroidesPatnodeLab_GTAC_VPL/BacteroidesPatnodePlate1_121323/LIB015096_22GMJHLT3_S14_L006_R1_001.fastq.gz | wc -l
  2861632 # 715,408 reads in this (by complete luck I randomly chose) B. theta sequencing data. Turns out I can get the # reads in plate1_info.csv, don't have to calculate this with zcat

# July 17, 2024
(biopython) [empchase@n0000 synthetic_contaminant]$ which python
/global/software/rocky-8.x86_64/python-3.11.6/python-3.11.6/python/3.11.6-gcc/bin/python ##points to system python executable, not python inside my environment
# created a test script to see if biopython works
(biopython) [empchase@n0000 synthetic_contaminant]$ python test_biopython.py 
  Biopython is not installed.
#install and reinstall
(biopython) [empchase@n0000 synthetic_contaminant]$ mamba remove biopython
(biopython) [empchase@n0000 synthetic_contaminant]$ mamba install -c conda-forge biopython
(biopython) [empchase@n0000 synthetic_contaminant]$ python test_biopython.py 
  Biopython is not installed.
#had to create a new environment for running biopython ("cmd_biopython") 
#because my "biopython" environment pointed to the savio python executable

# July 18, 2024
# making the synthetic contaminated file -- I've done a lot of work previously on these scripts but today it finally started working
# existing structure:
  # synthetic_contaminant/
    # two versions of contamination scripts. v2 is currently most up to date
    # test/ # contains pear test on ../subsetD4_read*.fastq.gz (as well as testing full length .fastq.gz files)
      # test_outputs #slurm and pear outputs and initial contaminate.py (renamed contaminate_v1.py) output
    # trimmomatic_test/ # contains script/dependent files/slurm outputs for trimmomatic
      cat /global/scratch/projects/fc_wolflab/bbioinfo/sarrah/batches/adapter_files/D4_adapters.fasta /global/scratch/projects/fc_wolflab/bbioinfo/sarrah/batches/adapter_files/F1_adapters.fasta > D4+F1_adapters.fasta
      # ^ adapter file for both adapters in concatenated fastq's
        #trimmomatic_test/ # trimmomatic outputs #i got lazy and named the subdirectory the same as the parent sorry
          mkdir fastqc
          # fastqc/ #fastqc outputs
            mkdir multiqc
            # multiqc/ #multiqc outputs
          mkdir fastqc_v0-11-9
          # fastqc_v0-11-9/ #corrected version to match srrah fastqc outputs
            mkdir multiqc
    mkdir spades_test
    # 
        


# July 19, 2024
# run spades on contaminated file

source activate bbio
# I couldn't figure out how to get spades to run but I found it in /global/scratch/projects/fc_wolflab/software so I ran
MODULEPATH=${MODULEPATH}:/global/scratch/projects/fc_wolflab/software
# this pulled up warnings when i ran module load xyz:
  Lmod Warning:  MODULEPATH directory: "/global/scratch/projects/fc_wolflab/software" has too many
  non-modulefiles (599). Please make sure that modulefiles are in their own directory and not mixed in with
  non-modulefiles (e.g. source code)
# to take the path out I ran this 
export MODULEPATH=$(echo $MODULEPATH | tr ':' '\n' | grep -v "/global/scratch/projects/fc_wolflab/software" | paste -sd:)
echo $MODULEPATH # no more thingy
  /global/software/rocky-8.x86_64/modfiles/langs:/global/software/rocky-8.x86_64/modfiles/tools:/global/software/rocky-8.x86_64/modfiles/compilers:/global/software/rocky-8.x86_64/modfiles/apps:/clusterfs/vector/home/groups/software/sl-7.x86_64/modfiles

MODULEPATH=${MODULEPATH}:/global/scratch/projects/fc_wolflab/software/modfiles

#create a new environment for spades, using python 3.8.8
conda create -n python3-8-8 python=3.8.8

# July 24, 2024
# test anvio pangenomics pipeline
# following along with anvio example
wget https://ndownloader.figshare.com/files/28834476 -O Prochlorococcus_31_genomes.tar.gz
tar -zxvf Prochlorococcus_31_genomes.tar.gz
cd Prochlorococcus_31_genomes
conda activate anvio-8
anvi-migrate *.db --migrate-dbs-safely
# keyboard interrupt bc I want the fast version
anvi-migrate *.db --migrate-quickly
  NEW MIGRATION TASK
  ===============================================
  Input file path ..............................: SS51.db
  Input file type ..............................: contigs
  Current Version ..............................: 20
  Target Version ...............................: 21
  Migration mode ...............................: Adventurous

  SQLite Version ...............................: 3.45.2

                                                                                                                                                                                                                    
  * Congratulations! Your contigs database is now version 21, which means it now
    contains two new empty tables. These tables are not as boring as they first
    appear, because you can now run `anvi-reaction-network` to store a network of
    the metabolic reactions that may be encoded by genes. A metabolic model
    representing the network can be exported from the database using `anvi-get-
    metabolic-model-file`.

#generate a genomes storage    
anvi-gen-genomes-storage -e external-genomes.txt -o PROCHLORO-GENOMES.db
  WARNING
  ===============================================
  Good news! Anvi'o found all these functions that are common to all of your
  genomes and will use them for downstream analyses and is very proud of you:
  'COG14_FUNCTION, COG14_CATEGORY'.

  Internal genomes .............................: 0 have been initialized.                                                                                                                                                         
  External genomes .............................: 31 found.                                                                                                                                                                        
                                                                                                                                                                                                                                  
  JUST FYI
  ===============================================
  Some of your genomes had gene calls identified by gene callers other than the
  anvi'o default, 'prodigal', and will not be processed. Use the `--debug` flag if
  this sounds important and you would like to see more of this message.

                                                                                                                                                                                                                                  
  * AS9601 is stored with 1,869 genes (0 of which were partial)
  * CCMP1375 is stored with 1,826 genes (0 of which were partial)                                                                                                                                                                  
  * EQPAC1 is stored with 1,892 genes (6 of which were partial)                                                                                                                                                                    
  * GP2 is stored with 1,825 genes (22 of which were partial)            
  The new genomes storage ......................: PROCHLORO-GENOMES.db (v7, signature: hash0cde9439)
  Number of genomes ............................: 31 (internal: 0, external: 31)
  Number of gene calls .........................: 60,223
  Number of partial gene calls .................: 409

# start pangenomics
touch pangenome.sh
sbatch pangenome.sh # runs the following code:
  anvi-pan-genome -g PROCHLORO-GENOMES.db \
                --project-name "Prochlorococcus_Pan" \
                --output-dir PROCHLORO \
                --num-threads 6 \
                --minbit 0.5 \
                --mcl-inflation 10 \
                --use-ncbi-blast
# 20 min run time
# scp everything to visualize interactively
# locally ran this to add misc data:
anvi-import-misc-data layer-additional-data.txt \
                      -p PROCHLORO/Prochlorococcus_Pan-PAN.db \
                      --target-data-table layers

  New layers additional data...
  ===============================================
  Data key "clade" .............................: Predicted type: str
  Data key "light" .............................: Predicted type: str

  New data added to the db for your layers .....: clade, light.


# July 25, 2024
# goal: get anvio running on some of our B theta genomes 
cd anvio_test
mkdir B_theta_test
#referenced the bbioinfo/plate1_info.csv to know which samples were b theta
touch btheta-genomes.txt #tab separated file wih b theta strain IDs_plate# and filepath


# August 2, 2024
# test anvio pipeline on vibrio tutorial

cd anvio_test
mkdir Vibrio_jascida
cd Vibrio_jascida

# download the pack
curl -L https://ndownloader.figshare.com/files/28965090 -o V_jascida_genomes.tar.gz
# unpack it
tar -zxvf V_jascida_genomes.tar.gz
# go into the new directory
cd V_jascida_genomes

# download fasta of reference genome (can filter for ref genome)
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002887615.1/
# use scp to move ncbi_dataset/data/GCF_002887615.1/GCF_002887615.1_ASM288761v1_genomic.fna into V_jascida_genomes
# from local terminal: scp /Users/emilychase/Desktop/WolfLab/Bacteroides_Bioinformatics/ncbi_dataset/ncbi_dataset/data/GCF_014131755.1/GCF_014131755.1_ASM1413175v1_genomic.fna empchase@dtn.brc.berkeley.edu:/global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/Vibrio_jascida/V_jascida_genomes
# OH NO I DOWNLOADED THE WRONG GENOME, not really sure what I just copied over, honestly it's probably a bacteroidetes genome
mv GCF_014131755.1_ASM1413175v1_genomic.fna ../..
# pull the correct genome
# local terminal: scp /Users/emilychase/Desktop/WolfLab/Bacteroides_Bioinformatics/anvio/V_jasicida_NCBI/ncbi_dataset/data/GCF_002887615.1/GCF_002887615.1_ASM288761v1_genomic.fna empchase@dtn.brc.berkeley.edu:/global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/Vibrio_jascida/V_jascida_genomes
mv GCF_002887615.1_ASM288761v1_genomic.fna ref_scaffolds.fasta

#for loop / prep
ls *fasta | awk 'BEGIN{FS="_"}{print $1}' > genomes.txt
for g in `cat genomes.txt`
do
    echo
    echo "Working on $g ..."
    echo
    anvi-script-reformat-fasta ${g}_scaffolds.fasta \
                               --min-len 2500 \
                               --simplify-names \
                               -o ${g}_scaffolds_2.5K.fasta
done
# if you want to look the command print statements, they're all on local anvio_log.sh

# generate contigs db
for g in `cat genomes.txt`
do
    echo
    echo "Working on $g ..."
    echo
    anvi-gen-contigs-database -f ${g}_scaffolds_2.5K.fasta \
                              -o V_jascida_${g}.db \
                              --num-threads 12 \
                              -n V_jascida_${g}
done

#Annotate contigs db

for g in *.db
do
    anvi-run-hmms -c $g --num-threads 18
    anvi-run-ncbi-cogs -c $g --num-threads 18 # error: need to run anvi-setup-ncbi-cogs first
    anvi-scan-trnas -c $g --num-threads 18
    anvi-run-scg-taxonomy -c $g --num-threads 18 # error: need to run anvi-setup-scg-taxonomy first
done

anvi-setup-ncbi-cogs --num-threads 18
anvi-setup-scg-taxonomy --num-threads 18

for g in *.db
do
    anvi-run-ncbi-cogs -c $g --num-threads 18
    anvi-run-scg-taxonomy -c $g --num-threads 18
done

# Tried this
anvi-display-contigs-stats *db
# realized it pulls up a window on the browser

# generate external genomes file -- automates creation of this tab-deliminated file from contigs dbs found in a directory
anvi-script-gen-genomes-file --input-dir . -o external-genomes.txt
  Contigs databases ............................: 8 found.                                                                                                                                        
  External genomes file ........................: external-genomes.txt

# investigate contamination
anvi-estimate-genome-completeness -e external-genomes.txt
  ValueError: numpy.dtype size changed, may indicate binary incompatibility. Expected 96 from C header, got 88 from PyObject

conda list numpy
  # packages in environment at /global/home/users/empchase/.conda/envs/anvio-8:
  #
  # Name                    Version                   Build  Channel
  numpy                     1.23.5                   pypi_0    pypi
conda list scikit-learn
  # packages in environment at /global/home/users/empchase/.conda/envs/anvio-8:
  #
  # Name                    Version                   Build  Channel
  scikit-learn              1.2.2                    pypi_0    pypi
# the instructions were to run these but I decided to set up new environment instead -- theoretically the numpy and scikitlearn from anvio should work together so maybe there was an anvio installation problem
# conda update numpy
# conda update scikit-learn

# try setting up new environment according to the folks online
conda --version                                                                                                                                          
  conda 23.3.1                     
conda activate conda
  EnvironmentNotWritableError: The current user does not have write permissions to the target environment.                                                                                        
  environment location: /global/software/sl-7.x86_64/modules/langs/python/3.10                                                                                                                  
  uid: 46620                                                                                                                                                                                    
  gid: 501 

conda create -y --name anvio-8-2 python=3.10

# August 7, 2024
# checkm installation
module load python
pip3 install checkm-genome
pip install hmmer


sbatch checkm_test.sh 
   Submitted batch job 20359262 # It seems that the CheckM data folder has not been set yet or has been removed. Please run 'checkm data setRoot'.


export PATH=$PATH:/global/home/users/empchase/.local/bin
checkm data setRoot

sbatch checkm_test.sh 
  Submitted batch job 20360183 #[2024-08-07 21:42:56] ERROR: Make sure prodigal is on your system path.

# edited checkm_test.sh to include module load prodigal
# export MODULEPATH=${MODULEPATH}:/clusterfs/vector/home/groups/software/sl-7.x86_64/modfiles
# module load prodigal
sbatch checkm_test.sh 
  Submitted batch job 20361193 #FileNotFoundError: [Errno 2] No such file or directory: '/global/home/users/empchase/.checkm/hmms/phylo.hmm'

# apparently I need to redo/go thru the checkm installation process myself -- not just pip install but some database installation steps too
curl -o ~/packages/checkm_data.tar.gz https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
    % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                  Dload  Upload   Total   Spent    Left  Speed
  100  275M  100  275M    0     0  6005k      0  0:00:46  0:00:46 --:--:-- 5641k
pushd ~/packages/
tar -xzvf checkm_data.tar.gz
# didn't extract into a directory oops
mkdir -p /global/scratch/projects/fc_wolflab/software/checkm_data
# move to wolf lab directory for lab use
mv distributions genome_tree hmms hmms_ssu img pfam selected_marker_sets.tsv taxon_marker_sets.tsv test_data .dmanifest /global/scratch/projects/fc_wolflab/software/checkm_data
export CHECKM_DATA_PATH=/global/scratch/projects/fc_wolflab/software/checkm_data
# alternatively could run checkm data setRoot /global/scratch/projects/fc_wolflab/software/checkm_data like I did before?

sbatch checkm_test.sh 
  Submitted batch job 20365698 # checkm is running now but: [2024-08-07 22:31:33] ERROR: Make sure pplacer is on your system path.
# added module load pplacer

sbatch checkm_test.sh 
  Submitted batch job 20366388

# August 8, 2024
# retry anvio dependency troubleshooting

cd bbioinfo/empchase/anvio_test/
conda activate anvio-8
pip show numpy
  Name: numpy
  Version: 1.23.5
  Summary: NumPy is the fundamental package for array computing with Python.
  Home-page: https://www.numpy.org
  Author: Travis E. Oliphant et al.
  Author-email: 
  License: BSD
  Location: /global/home/users/empchase/.conda/envs/anvio-8/lib/python3.10/site-packages
  Requires: 
  Required-by: anvio, biopython, contourpy, illumina-utils, matplotlib, numba, pandas, patsy, scikit-learn, scipy, seaborn, statsmodels

# iva said to change the version based on https://stackoverflow.com/questions/78650222/valueerror-numpy-dtype-size-changed-may-indicate-binary-incompatibility-expec 
pip install numpy==1.26.4
  Collecting numpy==1.26.4
    Downloading numpy-1.26.4-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (61 kB)
      ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 61.0/61.0 kB 1.9 MB/s eta 0:00:00
  Downloading numpy-1.26.4-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (18.2 MB)
    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 18.2/18.2 MB 35.1 MB/s eta 0:00:00
  Installing collected packages: numpy
    Attempting uninstall: numpy
      Found existing installation: numpy 1.23.5
      Uninstalling numpy-1.23.5:
        Successfully uninstalled numpy-1.23.5
  ERROR: pip's dependency resolver does not currently take into account all the packages that are installed. This behaviour is the source of the following dependency conflicts.
  anvio 8 requires numpy<=1.24, but you have numpy 1.26.4 which is incompatible.
  anvio 8 requires pandas==1.4.4, but you have pandas 2.2.2 which is incompatible.
  Successfully installed numpy-1.26.4'

# anvio 8 requires numpy<=1.24, but now I have numpy 1.26.4 which is incompatible. Going to retry with numpy 1.24.0
pip install numpy==1.24.0
  Collecting numpy==1.24.0
    Using cached numpy-1.24.0-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (5.6 kB)
  Downloading numpy-1.24.0-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (17.3 MB)
    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 17.3/17.3 MB 35.5 MB/s eta 0:00:00
  Installing collected packages: numpy
    Attempting uninstall: numpy
      Found existing installation: numpy 1.26.4
      Uninstalling numpy-1.26.4:
        Successfully uninstalled numpy-1.26.4
  ERROR: pip's dependency resolver does not currently take into account all the packages that are installed. This behaviour is the source of the following dependency conflicts.
  anvio 8 requires pandas==1.4.4, but you have pandas 2.2.2 which is incompatible.
  seaborn 0.13.2 requires numpy!=1.24.0,>=1.20, but you have numpy 1.24.0 which is incompatible.
  Successfully installed numpy-1.24.0'

conda env export > anvio-8-20240808-backup.yml

# August 9, 2024
conda activate anvio-8

pip install pandas==1.4.4
  Collecting pandas==1.4.4
    Using cached pandas-1.4.4-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (12 kB)
  Requirement already satisfied: python-dateutil>=2.8.1 in /global/home/users/empchase/.conda/envs/anvio-8/lib/python3.10/site-packages (from pandas==1.4.4) (2.9.0.post0)
  Requirement already satisfied: pytz>=2020.1 in /global/home/users/empchase/.conda/envs/anvio-8/lib/python3.10/site-packages (from pandas==1.4.4) (2024.1)
  Requirement already satisfied: numpy>=1.21.0 in /global/home/users/empchase/.conda/envs/anvio-8/lib/python3.10/site-packages (from pandas==1.4.4) (1.24.0)
  Requirement already satisfied: six>=1.5 in /global/home/users/empchase/.conda/envs/anvio-8/lib/python3.10/site-packages (from python-dateutil>=2.8.1->pandas==1.4.4) (1.16.0)
  Using cached pandas-1.4.4-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (11.6 MB)
  Installing collected packages: pandas
    Attempting uninstall: pandas
      Found existing installation: pandas 2.2.2
      Uninstalling pandas-2.2.2:
        Successfully uninstalled pandas-2.2.2
  ERROR: pip's dependency resolver does not currently take into account all the packages that are installed. This behaviour is the source of the following dependency conflicts.
  seaborn 0.13.2 requires numpy!=1.24.0,>=1.20, but you have numpy 1.24.0 which is incompatible.
  Successfully installed pandas-1.4.4'

conda remove numpy
pip uninstall numpy
conda install numpy=1.23.4

# retry the last anvio command
cd anvio_test/Vibrio_jascida/V_jascida_genomes
anvi-estimate-genome-completeness -e external-genomes.txt
  Traceback (most recent call last):
    File "/global/home/users/empchase/.conda/envs/anvio-8/bin/anvi-estimate-genome-completeness", line 7, in <module>
      import anvio
    File "/global/home/users/empchase/.conda/envs/anvio-8/lib/python3.10/site-packages/anvio/__init__.py", line 3562, in <module>
      from anvio.terminal import Run
    File "/global/home/users/empchase/.conda/envs/anvio-8/lib/python3.10/site-packages/anvio/terminal.py", line 12, in <module>
      import pandas as pd
    File "/global/home/users/empchase/.conda/envs/anvio-8/lib/python3.10/site-packages/pandas/__init__.py", line 22, in <module>
      from pandas.compat import is_numpy_dev as _is_numpy_dev
    File "/global/home/users/empchase/.conda/envs/anvio-8/lib/python3.10/site-packages/pandas/compat/__init__.py", line 15, in <module>
      from pandas.compat.numpy import (
    File "/global/home/users/empchase/.conda/envs/anvio-8/lib/python3.10/site-packages/pandas/compat/numpy/__init__.py", line 4, in <module>
      from pandas.util.version import Version
    File "/global/home/users/empchase/.conda/envs/anvio-8/lib/python3.10/site-packages/pandas/util/__init__.py", line 1, in <module>
      from pandas.util._decorators import (  # noqa:F401
    File "/global/home/users/empchase/.conda/envs/anvio-8/lib/python3.10/site-packages/pandas/util/_decorators.py", line 14, in <module>
      from pandas._libs.properties import cache_readonly  # noqa:F401
    File "/global/home/users/empchase/.conda/envs/anvio-8/lib/python3.10/site-packages/pandas/_libs/__init__.py", line 13, in <module>
      from pandas._libs.interval import Interval
    File "pandas/_libs/interval.pyx", line 1, in init pandas._libs.interval
    File "pandas/_libs/hashtable.pyx", line 1, in init pandas._libs.hashtable
    File "pandas/_libs/missing.pyx", line 1, in init pandas._libs.missing
    File "/global/home/users/empchase/.conda/envs/anvio-8/lib/python3.10/site-packages/pandas/_libs/tslibs/__init__.py", line 30, in <module>
      from pandas._libs.tslibs.conversion import (
    File "pandas/_libs/tslibs/conversion.pyx", line 1, in init pandas._libs.tslibs.conversion
    File "pandas/_libs/tslibs/timezones.pyx", line 14, in init pandas._libs.tslibs.timezones
    File "/global/home/users/empchase/.conda/envs/anvio-8/lib/python3.10/site-packages/dateutil/tz/__init__.py", line 2, in <module>
      from .tz import *
    File "/global/home/users/empchase/.conda/envs/anvio-8/lib/python3.10/site-packages/dateutil/tz/tz.py", line 19, in <module>
      import six
  ModuleNotFoundError: No module named 'six'

conda install --channel "vpython" six
  Collecting package metadata (current_repodata.json): done

anvi-estimate-genome-completeness -e external-genomes.txt
  Traceback (most recent call last):
    File "/global/home/users/empchase/.conda/envs/anvio-8/bin/anvi-estimate-genome-completeness", line 8, in <module>
      import anvio.utils as utils
    File "/global/home/users/empchase/.conda/envs/anvio-8/lib/python3.10/site-packages/anvio/utils.py", line 28, in <module>
      import Bio.PDB as PDB
  ModuleNotFoundError: No module named 'Bio'
# export MODULEPATH=${MODULEPATH}:/clusterfs/vector/home/groups/software/sl-7.x86_64/modfiles didn't resolve the issue

mamba install -c conda-forge biopython

anvi-estimate-genome-completeness -e external-genomes.txt
  /global/home/users/empchase/.conda/envs/anvio-8/lib/python3.10/site-packages/numba/__init__.py:48: UserWarning: A NumPy version >=1.23.5 and <2.3.0 is required for this version of SciPy (detected version 1.23.4)
    import scipy
  +---------------+----------+--------------+----------------+----------------+--------------+----------------+
  | genome name   | domain   |   confidence |   % completion |   % redundancy |   num_splits |   total length |
  +===============+==========+==============+================+================+==============+================+
  | V_jascida_12  | BACTERIA |            1 |            100 |           5.63 |          282 |        5814827 |
  +---------------+----------+--------------+----------------+----------------+--------------+----------------+
  | V_jascida_13  | BACTERIA |            1 |            100 |           5.63 |          283 |        5812569 |
  +---------------+----------+--------------+----------------+----------------+--------------+----------------+
  | V_jascida_14  | BACTERIA |            1 |            100 |           5.63 |          290 |        5861661 |
  +---------------+----------+--------------+----------------+----------------+--------------+----------------+
  | V_jascida_47  | BACTERIA |            1 |            100 |           5.63 |          282 |        5821176 |
  +---------------+----------+--------------+----------------+----------------+--------------+----------------+
  | V_jascida_52  | BACTERIA |            1 |            100 |          21.13 |          362 |        6303185 |
  +---------------+----------+--------------+----------------+----------------+--------------+----------------+
  | V_jascida_53  | BACTERIA |            1 |            100 |           5.63 |          289 |        5853972 |
  +---------------+----------+--------------+----------------+----------------+--------------+----------------+
  | V_jascida_55  | BACTERIA |            1 |            100 |           5.63 |          289 |        5857980 |
  +---------------+----------+--------------+----------------+----------------+--------------+----------------+
  | V_jascida_ref | BACTERIA |            1 |            100 |           5.63 |          298 |        5989523 |
  +---------------+----------+--------------+----------------+----------------+--------------+----------------+

  * The 'domain' shown in this table is the domain anvi'o predicted for your contigs
    in a given bin with the amount of confidence for that call in the 'domain'
    column. If the domain is 'mixed', it means it is very likely you have contigs
    in your bin that spans accross multiple domains of life (maybe it is worth a
    Nobel, but more likely it is just garbage). If the domain is 'blank', it means
    anvi'o did not find enough signal from single-copy core genes to determine
    what domain it should use to estimate the completion and redundancy of this
    bin accurately. You can get much more information about these estimations by
    running the same command line with the additinal flag `--debug`.

# pause on this environment, try new one
conda deactivate
conda activate anvio-8-2
anvi-estimate-genome-completeness -e external-genomes.txt

anvi-profile -c V_jascida_52.db \
             --sample-name V_jascida_52 \
             --output-dir V_jascida_52 \
             --blank

  input_bam .........................................: None                               
  output_dir ........................................: /global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/Vibrio_jascida/V_jascida_genomes/V_jascida_52
  num_contigs_after_M ...............................: 201
  num_splits ........................................: 362
  total_length ......................................: 6,303,185
                                                                                          
  Additional data added to the new profile DB .......: total_reads_mapped
  New items order ...................................: "tnf-splits:euclidean:ward" (type newick) has been added to the database...
  New items order ...................................: "tnf:euclidean:ward" (type newick) has been added to the database...

  * Happy 😇


  ✓ anvi-profile took 0:00:25.499911

# on local, ran: 
#scp empchase@dtn.brc.berkeley.edu:/global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/Vibrio_jascida/V_jascida_genomes/V_jascida_52.db .
# didn't work

# does it work on savio?
anvi-interactive -c V_jascida_52.db                  -p V_jascida_52/PROFILE.db
Contigs DB ...................................: Initialized: V_jascida_52.db (v. 21)
Interactive mode .............................: full                        
Auxiliary Data ...............................: Found: V_jascida_52/AUXILIARY-DATA.db (v. 2)
Profile Super ................................: Initialized with all 362 splits: V_jascida_52/PROFILE.db (v. 38)
                                                                            
WARNING
===============================================
Even though the following layer names were associated with your data, they will
not be shown in your interactive display since none of the splits you are in
terested at this point had any hits in those layers: "Ribosomal_RNA_12S,
Ribosomal_RNA_18S, Ribosomal_RNA_28S, Ribosomal_RNA_5S". If you would like every
layer to be shown in the interface even when there are nohits, please include
the flag `--show-all-layers` in your command.

Contigs DB ...................................: V_jascida_52.db
Metagenome mode ..............................: False
                                                                            
* The server is up and running 🎉

WARNING
===============================================
If you are using OSX and if the server terminates prematurely before you can see
anything in your browser, try running the same command by putting 'sudo ' at the
beginning of it (you will be prompted to enter your password if sudo requires
super user credentials on your system). If your browser does not show up, try
manually entering the URL shown below into the address bar of your favorite
browser. *cough* CHROME *cough*.


Server address ...............................: http://0.0.0.0:8080

* When you are ready, press CTRL+C once to terminate the server and go back to the
  command line.

^C
* The server is being terminated...
# when savio launches the url it just says "Not Found."


# oh on local I had forgotten the complementary databases
#scp -r empchase@dtn.brc.berkeley.edu:/global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/Vibrio_jascida/V_jascida_genomes/V_jascida_52 .
anvi-interactive -c V_jascida_52.db                  -p V_jascida_52/PROFILE.db
# launches the appropriate page/interface. 
# Draw image, select branches of the tree based on what the contamination looks like into two bins (contam, clean),
# click "store bin collection" 
# I called my bin collection default_20240809_EClocal

#Bring the updated db's back -- first need to move old db's
mkdir localclean
mv V_jascida_52* localclean/
# scp V_jascida_52.db empchase@dtn.brc.berkeley.edu:/global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/Vibrio_jascida/V_jascida_genomes
# scp -r V_jascida_52 empchase@dtn.brc.berkeley.edu:/global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/Vibrio_jascida/V_jascida_genomes

anvi-split -p V_jascida_52/PROFILE.db -c V_jascida_52.db -C default_20240809_EClocal -o V_jascida_52_SPLIT

sed 's/V_jascida_52.db/V_jascida_52_SPLIT\/V_jascida_52_CLEAN\/CONTIGS.db/g' external-genomes.txt > external-genomes-final.txt

anvi-gen-genomes-storage -e external-genomes-final.txt \
                         -o V_jascida-GENOMES.db

#normal outputs
anvi-pan-genome -g V_jascida-GENOMES.db \
                --project-name V_jascida \
                --num-threads 8

#slurm died in the middle but this was good progress!
  [09 Aug 24 16:57:22 Aligning amino acid sequences for genes in gene (...) ETA: 14sslurmstepd: error: *** STEP 20404816.0 ON n0030.savio2 CANCELLED AT 2024-08-09T16:57:38 ***
  [09 Aug 24 16:57:22 Aligning amino acid sequences for genes in gene (...) ETA: 14ssrun: Job step aborted: Waiting up to 32 seconds for job step to finish.
  [09 Aug 24 16:58:07 Vectors from 2 matrices] Clustering data with "ward" lin (...)srun: error: n0030.savio2: task 0: Killed
  srun: launch/slurm: _step_signal: Terminating StepId=20404816.0
  srun: error: Timed out waiting for job step to complete

# August 12, 2024
cd bbioinfo/empchase/anvio_test/Vibrio_jascida/V_jascida_genomes

anvi-pan-genome -g V_jascida-GENOMES.db \
                --project-name V_jascida \
                --num-threads 12
# not running in interactive for some reason?
^C
# try submitting as a job?
touch pan-genome_test.sh
sbatch pan-genome_test.sh 
  Submitted batch job 20405856

#checkm
# in interactive mode
cd /global/scratch/projects/fc_wolflab/bbioinfo/empchase/checkm_test
bash checkm_test.sh
  /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/pplacer/1.1b19:/clusterfs/vector/home/groups/software/sl-7.x86_64/modules/prodigal/2.6.3:/global/home/users/empchase/bin:/global/software/rocky-8.x86_64/manual/modules/tools/sq/0.1.0/bin:/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-8.5.0/emacs-29.1-d5gecqgjml5m6y3gny5sqwqxy4mys4ex/bin:/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-8.5.0/code-server-4.91.1-rywu6um26xekdmcnat7o4xzuwfzokfls/lib/vscode/bin/remote-cli:/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-8.5.0/code-server-4.91.1-rywu6um26xekdmcnat7o4xzuwfzokfls/bin:/global/software/sl-7.x86_64/modules/langs/python/3.10/condabin:/usr/share/Modules/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/global/home/groups/allhands/bin:/global/home/users/empchase/.local/bin
  [2024-08-12 18:17:56] INFO: CheckM v1.2.3
  [2024-08-12 18:17:56] INFO: checkm lineage_wf -x fasta /global/scratch/projects/fc_wolflab/bbioinfo/sarrah/batches/trimmed/Plate1/A1/100/assembly /global/scratch/projects/fc_wolflab/bbioinfo/empchase/checkm_test/checkm_results_rocky8
  [2024-08-12 18:17:56] INFO: CheckM data: /global/scratch/projects/fc_wolflab/software/checkm_data
  [2024-08-12 18:17:56] INFO: [CheckM - tree] Placing bins in reference genome tree.
  [2024-08-12 18:17:56] INFO: Identifying marker genes in 3 bins with 1 threads:
      Finished processing 3 of 3 (100.00%) bins.
  [2024-08-12 18:20:41] INFO: Saving HMM info to file.
  [2024-08-12 18:20:41] INFO: Calculating genome statistics for 3 bins with 1 threads:
      Finished processing 3 of 3 (100.00%) bins.
  [2024-08-12 18:20:42] INFO: Extracting marker genes to align.
  [2024-08-12 18:20:42] INFO: Parsing HMM hits to marker genes:
      Finished parsing hits for 3 of 3 (100.00%) bins.
  [2024-08-12 18:20:42] INFO: Extracting 43 HMMs with 1 threads:
      Finished extracting 43 of 43 (100.00%) HMMs.
  [2024-08-12 18:20:45] INFO: Aligning 43 marker genes with 1 threads:
      Finished aligning 43 of 43 (100.00%) marker genes.
  [2024-08-12 18:20:53] ERROR: Make sure pplacer is on your system path.

# same error as August 7
which pplacer
  /usr/bin/which: no pplacer in (/clusterfs/vector/home/groups/software/sl-7.x86_64/modules/pplacer/1.1b19:/global/home/users/empchase/bin:/global/software/rocky-8.x86_64/manual/modules/tools/sq/0.1.0/bin:/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-8.5.0/emacs-29.1-d5gecqgjml5m6y3gny5sqwqxy4mys4ex/bin:/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-8.5.0/code-server-4.91.1-rywu6um26xekdmcnat7o4xzuwfzokfls/lib/vscode/bin/remote-cli:/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-8.5.0/code-server-4.91.1-rywu6um26xekdmcnat7o4xzuwfzokfls/bin:/global/software/sl-7.x86_64/modules/langs/python/3.10/condabin:/usr/share/Modules/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/global/home/groups/allhands/bin)

module show pplacer
  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    /clusterfs/vector/home/groups/software/sl-7.x86_64/modfiles/pplacer/1.1a19:
  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  whatis("loads the environment for pplacer 1.1b19")
  setenv("PPLACER_DIR","/clusterfs/vector/home/groups/software/sl-7.x86_64/modules/pplacer/1.1b19")
  prepend_path("PATH","/clusterfs/vector/home/groups/software/sl-7.x86_64/modules/pplacer/1.1b19")
  help([[loads the environment for pplacer 1.1b19
  ]])

ls -l /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/pplacer/1.1b19
ls: cannot access '/clusterfs/vector/home/groups/software/sl-7.x86_64/modules/pplacer/1.1b19': No such file or directory
[empchase@n0130 checkm_test]$ ls -l /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/pplacer/
total 4
drwxrwxr-x 1 guohaoang cgrlsw 196 Nov 20  2017 1.1a19

# change checkm_test.sh
# module load pplacer/1.1a19

bash checkm_test.sh
# same error: ERROR: Make sure pplacer is on your system path.

# try installing pplacer on my own:
conda install bioconda::pplacer
# permission denied



# August 15, 2024
# checkm
# amogh fixed pplacer !!
cd /global/scratch/projects/fc_wolflab/bbioinfo/empchase/checkm_test
bash checkm_test.sh 
  /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/pplacer/1.1a19:/clusterfs/vector/home/groups/software/sl-7.x86_64/modules/prodigal/2.6.3:/global/home/users/empchase/bin:/global/software/rocky-8.x86_64/manual/modules/tools/sq/0.1.0/bin:/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-8.5.0/emacs-29.1-d5gecqgjml5m6y3gny5sqwqxy4mys4ex/bin:/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-8.5.0/code-server-4.91.1-rywu6um26xekdmcnat7o4xzuwfzokfls/lib/vscode/bin/remote-cli:/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-8.5.0/code-server-4.91.1-rywu6um26xekdmcnat7o4xzuwfzokfls/bin:/global/software/sl-7.x86_64/modules/langs/python/3.10/condabin:/usr/share/Modules/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/global/home/groups/allhands/bin:/global/home/users/empchase/.local/bin
  -------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Bin Id            Marker lineage         # genomes   # markers   # marker sets   0    1    2   3   4   5+   Completeness   Contamination   Strain heterogeneity  
  -------------------------------------------------------------------------------------------------------------------------------------------------------------------
    scaffolds   o__Bacteroidales (UID2657)      160         491           269        2   481   7   1   0   0       99.26            2.01              30.00          
    contigs     o__Bacteroidales (UID2657)      160         491           269        2   481   7   1   0   0       99.26            2.01              30.00          
    before_rr   o__Bacteroidales (UID2657)      160         491           269        2   481   7   1   0   0       99.26            2.01              30.00          
  -------------------------------------------------------------------------------------------------------------------------------------------------------------------
  [2024-08-15 14:58:22] INFO: { Current stage: 0:00:02.421 || Total: 0:20:39.589 }

# edited the checkm_contaminant.sh script to reflect what worked for checkm_test
# didn't put in print path statement
sbatch checkm_contaminant.sh 
  Submitted batch job 20457544


# back to anvio
cd /global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/Vibrio_jascida/V_jascida_genomes
conda activate anvio-8-2

anvi-pan-genome -g V_jascida-GENOMES.db \
                --project-name V_jascida \
                --num-threads 8

#now it's working again, mustve been some system issue
File/Path Error: The output file '/global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_tes
                 t/Vibrio_jascida/V_jascida_genomes/V_jascida/V_jascida-PAN.db' already exists.  
                 Generally speaking anvi'o tries to avoid overwriting stuff. But you can always  
                 use the flag `--force-overwrite` to instruct anvi'o to delete the existing file 
                 first.  

anvi-pan-genome -g V_jascida-GENOMES.db \
    --project-name V_jascida \
    --num-threads 8 --force-overwrite


anvi-compute-genome-similarity --external-genomes external-genomes-final.txt \
                               --program pyANI \ 
                               --output-dir ANI \
                               --num-threads 6 \
                               --pan-db V_jascida/V_jascida-PAN.db 
  CITATION
  ===============================================
  Anvi'o will use 'PyANI' by Pritchard et al. (DOI: 10.1039/C5AY02550H) to compute
  ANI. If you publish your findings, please do not forget to properly credit their
  work.

  [PyANI] Num threads to use ...................: 6
  [PyANI] Alignment method .....................: ANIb
  [PyANI] Log file path ........................: /tmp/tmp4479gwwd

                                                                                                                                                                                                                                      
  ANI RESULTS
  ===============================================
  * Matrix and clustering of 'alignment coverage' written to output directory
  * Matrix and clustering of 'alignment lengths' written to output directory
  * Matrix and clustering of 'hadamard' written to output directory
  * Matrix and clustering of 'percentage identity' written to output directory
  * Matrix and clustering of 'similarity errors' written to output directory
  * Matrix and clustering of 'full percentage identity' written to output directory

  MISC DATA MAGIC FOR YOUR PAN DB
  ===============================================
  * Additional data and order for ANI alignment coverage are now in pan db                                                                                                                                                             
  * Additional data and order for ANI alignment lengths are now in pan db                                                                                                                                                              
  * Additional data and order for ANI hadamard are now in pan db                                                                                                                                                                       
  * Additional data and order for ANI percentage identity are now in pan db                                                                                                                                                            
  * Additional data and order for ANI similarity errors are now in pan db                                                                                                                                                              
  * Additional data and order for ANI full percentage identity are now in pan db                                                                                                                                                       

  ✓ anvi-compute-genome-similarity took 0:06:28.552056    '                         

# on local
  # # running savio generated pan genome (+ANI)

  # anvi-display-pan -p V_jascida-PAN.db \
  #                 -g V_jascida-GENOMES.db# 


# August 16, 2024
# try out on b theta
# previously 
# anvi-script-reformat-fasta $contig_path -o ${sample}_assembly_contigs-fasta --simplify-names --prefix $sample   
# anvi-gen-contigs-database -f ${sample}_assembly_contigs-fasta -o ../${sample}_contigs.db --external-gene-calls ${sample}_gene_calls.txt --project-name ${sample}_ContigDB --ignore-internal-stop-codons
cd /global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/B_theta_test

anvi-setup-ncbi-cogs
anvi-setup-scg-taxonomy --num-threads 18

for g in *.db
do
    anvi-run-hmms -c $g --num-threads 18
    anvi-run-ncbi-cogs -c $g --num-threads 18 
    anvi-scan-trnas -c $g --num-threads 18
    anvi-run-scg-taxonomy -c $g --num-threads 18
done

anvi-profile -c A4_75_WH502_contigs.db --sample-name A8_VPI_C11_15 --output-dir A8_VPI_C11_15 --blank
# wrong labels
rm -r A8_VPI_C11_15
anvi-profile -c A4_75_WH502_contigs.db --sample-name A4_75_WH502 --output-dir A4_75_WH502 --blank
anvi-profile -c A8_VPI_C11_15_contigs.db --sample-name A8_VPI_C11_15 --output-dir A8_VPI_C11_15 --blank
anvi-profile -c C4_WH507_contigs.db --sample-name C4_WH507 --output-dir C4_WH507 --blank

mv betheta_prep ..
# run local
#scp -r empchase@dtn.brc.berkeley.edu:/global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/B_theta_test .
# anvi-interactive -c A4_75_WH502_contigs.db -p A4_75_WH502/PROFILE.db
# need to split ^
# anvi-interactive -c A8_VPI_C11_15_contigs.db -p A8_VPI_C11_15/PROFILE.db
# completeness estimate is 97% when i take out shortest splits, 98% when I leave them in. 
# anvi-interactive -c C4_WH507_contigs.db -p C4_WH507/PROFILE.db
# need to split ^
mkdir not_cleaned
(anvio-8-2) [empchase@n0204 B_theta_test]$ mv A4* not_cleaned/
(anvio-8-2) [empchase@n0204 B_theta_test]$ mv C4_WH507* not_cleaned/
# local: scp -r A4_75_WH502* empchase@dtn.brc.berkeley.edu:/global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/B_theta_test
# scp -r C4_WH507* empchase@dtn.brc.berkeley.edu:/global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/B_theta_test

anvi-split -p A4_75_WH502/PROFILE.db \
           -c A4_75_WH502_contigs.db \
           -C default \
           -o A4_75_WH502_SPLIT

anvi-split -p C4_WH507/PROFILE.db \
           -c C4_WH507_contigs.db \
           -C default \
           -o C4_WH507_SPLIT

sed -e 's/A4_75_WH502_contigs.db/A4_75_WH502_SPLIT\/Clean\/CONTIGS.db/g' -e 's/C4_WH507_contigs.db/C4_WH507_SPLIT\/Clean\/CONTIGS.db/g' external-genomes.txt > external-genomes-final.txt

# next to run
anvi-gen-genomes-storage -e external-genomes-final.txt \
                         -o Btheta_test-GENOMES.db
                         JUST FYI
  ===============================================
  Some of your genomes had gene calls identified by gene callers other than the
  anvi'o default, 'prodigal', and will not be processed. Use the `--debug` flag if
  this sounds important and you would like to see more of this message.

                                                                                                                                                                                                                        

  Config Error: None of your genomes seem to have a gene call, which is a typical error you get
                if you are working with contigs databases with external gene calls. You can    
                solve it by looking at the output of the program `anvi-db-info` for a given    
                contigs database in your collection, and use one of the gene caller sources    
                listed in the output using the `--gene-caller` parameter. '

anvi-db-info A8_VPI_C11_15_contigs.db 

  DB Info (no touch)
  ===============================================
  Database Path ................................: A8_VPI_C11_15_contigs.db
  description ..................................: [Not found, but it's OK]'
  db_type ......................................: contigs (variant: unknown)
  version ......................................: 21


  DB Info (no touch also)
  ===============================================
  project_name .................................: A8_VPI_C11_15_ContigDB
  contigs_db_hash ..............................: hash38b50b68
  split_length .................................: 20000
  kmer_size ....................................: 4
  num_contigs ..................................: 660
  total_length .................................: 6590518
  num_splits ...................................: 895
  gene_level_taxonomy_source ...................: None
  genes_are_called .............................: 1
  external_gene_calls ..........................: 1
  external_gene_amino_acid_seqs ................: 0
  skip_predict_frame ...........................: 0
  splits_consider_gene_calls ...................: 1
  trna_taxonomy_was_run ........................: 0
  trna_taxonomy_database_version ...............: None
  creation_date ................................: 1721962969.24523
  gene_function_sources ........................: COG20_PATHWAY,Transfer_RNAs,COG20_FUNCTION,COG20_CATEGORY
  scg_taxonomy_was_run .........................: 1
  scg_taxonomy_database_version ................: v214.1

  * Please remember that it is never a good idea to change these values. But in some
    cases it may be absolutely necessary to update something here, and a
    programmer may ask you to run this program and do it. But even then, you
    should be extremely careful.


  AVAILABLE GENE CALLERS
  ===============================================
  * 'Transfer_RNAs' (71 gene calls)
  * 'Ribosomal_RNA_23S' (1 gene calls)
  * 'Ribosomal_RNA_16S' (1 gene calls)
  * 'Prodigal' (5,044 gene calls)


  AVAILABLE FUNCTIONAL ANNOTATION SOURCES
  ===============================================
  * COG20_CATEGORY (3,402 annotations)
  * COG20_FUNCTION (3,402 annotations)
  * COG20_PATHWAY (711 annotations)
  * Transfer_RNAs (71 annotations)


  AVAILABLE HMM SOURCES
  ===============================================
  * 'Archaea_76' (76 models with 38 hits)
  * 'Bacteria_71' (71 models with 73 hits)
  * 'Protista_83' (83 models with 2 hits)
  * 'Ribosomal_RNA_12S' (1 model with 0 hits)
  * 'Ribosomal_RNA_16S' (3 models with 1 hit)
  * 'Ribosomal_RNA_18S' (1 model with 0 hits)
  * 'Ribosomal_RNA_23S' (2 models with 1 hit)
  * 'Ribosomal_RNA_28S' (1 model with 0 hits)
  * 'Ribosomal_RNA_5S' (5 models with 0 hits)
  * 'Transfer_RNAs' (61 models with 71 hits)

# so there are gene calls ? Maybe try one of the other sources? On debug mode?
anvi-gen-genomes-storage -e external-genomes-final.txt -o Btheta_test-GENOMES.db --gene-caller COG20_CATEGORY --debug
  PLEASE READ CAREFULLY
  ===============================================
  Some of your genomes had gene calls identified by gene callers other than the
  gene caller anvi'o used, which was set to 'COG20_CATEGORY' either by default, or
  because you asked for it. The following genomes contained genes that were not
  processed (this may be exactly what you expect to happen, but if was not, you
  may need to use the `--gene-caller` flag to make sure anvi'o is using the gene
  caller it should be using): A8_VPI_C11_15_ContigDB (5044 gene calls by
  "Prodigal", 1 gene calls by "Ribosomal_RNA_23S", 1 gene calls by
  "Ribosomal_RNA_16S", 71 gene calls by "Transfer_RNAs"), C4_WH507_ContigDB (4716
  gene calls by "Prodigal", 62 gene calls by "Transfer_RNAs"), btheta_WH502_test
  (4713 gene calls by "Prodigal", 1 gene calls by "Ribosomal_RNA_23S", 1 gene
  calls by "Ribosomal_RNA_16S", 65 gene calls by "Transfer_RNAs").
  Config Error: None of your genomes seem to have a gene call, which is a typical error you get
              if you are working with contigs databases with external gene calls. You can    
              solve it by looking at the output of the program `anvi-db-info` for a given    
              contigs database in your collection, and use one of the gene caller sources    
              listed in the output using the `--gene-caller` parameter. 

# Try again with prodigal, which I thought was the default but idk anything anymore
anvi-gen-genomes-storage -e external-genomes-final.txt -o Btheta_test-GENOMES.db --gene-caller Prodigal --debug
  PLEASE READ CAREFULLY
  ===============================================
  Some of your genomes had gene calls identified by gene callers other than the
  gene caller anvi'o used, which was set to 'Prodigal' either by default, or
  because you asked for it. The following genomes contained genes that were not
  processed (this may be exactly what you expect to happen, but if was not, you
  may need to use the `--gene-caller` flag to make sure anvi'o is using the gene
  caller it should be using): A8_VPI_C11_15_ContigDB (1 gene calls by
  "Ribosomal_RNA_23S", 1 gene calls by "Ribosomal_RNA_16S", 71 gene calls by
  "Transfer_RNAs"), C4_WH507_ContigDB (62 gene calls by "Transfer_RNAs"),
  btheta_WH502_test (1 gene calls by "Ribosomal_RNA_23S", 1 gene calls by
  "Ribosomal_RNA_16S", 65 gene calls by "Transfer_RNAs").

                                                                                                                                                                                                                        
  * A8_VPI_C11_15_ContigDB is stored with 5,044 genes (0 of which were partial)
  * C4_WH507_ContigDB is stored with 4,716 genes (0 of which were partial)                                                                                                                                               
  * btheta_WH502_test is stored with 4,713 genes (0 of which were partial)                                                                                                                                               

  The new genomes storage ......................: Btheta_test-GENOMES.db (v7, signature: hash6de31a02)
  Number of genomes ............................: 3 (internal: 0, external: 3)
  Number of gene calls .........................: 14,473
  Number of partial gene calls .................: 0
# okay it worked!
# then compute pangenome
anvi-pan-genome -g Btheta_test-GENOMES.db \
                --project-name Btheta_test \
                --num-threads 12

# calculate ANI to annotate the pangenome with! Uses pyANI
anvi-compute-genome-similarity --external-genomes external-genomes-final.txt --program pyANI --output-dir ANI --num-threads 18 --pan-db Btheta_test/Btheta_test-PAN.db
# i think it was okay but i lost the output
# also just realized I forgot to put the reference genome along with everyone else ...
# we'll start here though and put the ref genome in momentarily https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_014131755.1/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED

#for now scp pangenome over
# scp empchase@dtn.brc.berkeley.edu:/global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/B_theta_test/Btheta_test-GENOMES.db .
# scp empchase@dtn.brc.berkeley.edu:/global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/B_theta_test/Btheta_test/Btheta_test-PAN.db .
# anvi-display-pan -p Btheta_test-PAN.db -g Btheta_test-GENOMES.db

# it worked! It is a little meaningless without the reference genome, so let's do this workflow one more time from the beginning alongside the ref genome
cd /global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/
mv btheta_prep B_theta_test
mkdir B_theta_test_2
cd B_theta_test_2

# local scp ncbi_dataset/ncbi_dataset/data/GCA_014131755.1/GCA_014131755.1_ASM1413175v1_genomic.fna empchase@dtn.brc.berkeley.edu:/global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/B_theta_test_2
# curious what this does
curl -L https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_014131755.1/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED -o bthetancbi_genome_test.fna.zip --compressed
  Warning: Binary output can mess up your terminal. Use "--output -" to tell 
  Warning: curl to output it to your terminal anyway, or consider "--output 
  Warning: <FILE>" to save to a file.
# I'm afraid of what they mean by mess up your terminal so I'm just gonna stop trying there
cp /global/scratch/projects/fc_wolflab/bbioinfo/sarrah/batches/trimmed_all/Plate1/C4/B_theta_WH507_assembly/contigs.fasta .
cp /global/scratch/projects/fc_wolflab/bbioinfo/sarrah/batches/trimmed_all/Plate1/A4/B_theta_WH502_assembly/contigs.fasta .
cp /global/scratch/projects/fc_wolflab/bbioinfo/sarrah/batches/trimmed_all/Plate1/A8/B_theta_VPI-C11-15_assembly/contigs.fasta .
# these commands just rename "contigs.fasta" over and over lol
cp /global/scratch/projects/fc_wolflab/bbioinfo/sarrah/batches/trimmed_all/Plate1/C4/B_theta_WH507_assembly/contigs.fasta ./B_theta_WH507-contigs.fasta
cp /global/scratch/projects/fc_wolflab/bbioinfo/sarrah/batches/trimmed_all/Plate1/A4/B_theta_WH502_assembly/contigs.fasta ./B_theta_WH502-contigs.fasta
cp /global/scratch/projects/fc_wolflab/bbioinfo/sarrah/batches/trimmed_all/Plate1/A8/B_theta_VPI-C11-15_assembly/contigs.fasta ./B_theta_VPI-C11-15-contigs.fasta


# checkm on reference genome
cd /global/scratch/projects/fc_wolflab/bbioinfo/empchase/checkm_test
touch checkm_reference.sh # directed to the fna file I just got from ncbi
sbatch checkm_reference.sh 
  Submitted batch job 20491602 #failed: i think bc I put the wrong extension on the input (fasta instead of fna)
sbatch checkm_reference.sh 
  Submitted batch job 20491606

# August 20, 2024
# try out taxonomic assignment on B theta genomes
cd /global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/B_theta_test_2
conda activate anvio-8-2

mv GCA_014131755.1_ASM1413175v1_genomic.fna ref-contigs.fasta # b theta reference genome https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_014131755.1/ 

for g in *.fasta
> do
> echo
> echo "Working on $g ..."
> echo
> anvi-gen-contigs-database -f $g -o B_theta_${g}.db --num-threads 18 -n B_theta_${g}_name
> done
    Working on B_theta_VPI-C11-15-contigs.fasta ...

    Input FASTA file .............................: /global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/B_theta_test_2/B_theta_VPI-C11-15-contigs.fasta
    Name .........................................: B_theta_B_theta_VPI-C11-15-contigs.fasta_name
    Description ..................................: No description is given
    Num threads for gene calling .................: 18                 

    Finding ORFs in contigs
    ===============================================
    Genes ........................................: /tmp/tmptkkmq5yk/contigs.genes
    Amino acid sequences .........................: /tmp/tmptkkmq5yk/contigs.amino_acid_sequences
    Log file .....................................: /tmp/tmptkkmq5yk/00_log.txt

    CITATION
    ===============================================
    Anvi'o will use 'prodigal' by Hyatt et al (doi:10.1186/1471-2105-11-119) to
    identify open reading frames in your data. When you publish your findings,
    please do not forget to properly credit their work.

    Result .......................................: Prodigal (v2.6.3) has identified 5665 genes.

                                                                      
    CONTIGS DB CREATE REPORT
    ===============================================
    Split Length .................................: 20,000
    K-mer size ...................................: 4
    Skip gene calling? ...........................: False
    External gene calls provided? ................: False
    Ignoring internal stop codons? ...............: False
    Splitting pays attention to gene calls? ......: True
    Contigs with at least one gene call ..........: 642 of 660 (97.3%) 
    Contigs database .............................: A new database, B_theta_B_theta_VPI-C11-15-contigs.fasta.db, has been created.
    Number of contigs ............................: 660
    Number of splits .............................: 895
    Total number of nucleotides ..................: 6,590,518
    Gene calling step skipped ....................: False
    Splits broke genes (non-mindful mode) ........: False
    Desired split length (what the user wanted) ..: 20,000
    Average split length (what anvi'o gave back) .: 21,401

    ✓ anvi-gen-contigs-database took 0:00:16.428835

    Working on B_theta_WH502-contigs.fasta ...

    Input FASTA file .............................: /global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/B_theta_test_2/B_theta_WH502-contigs.fasta
    Name .........................................: B_theta_B_theta_WH502-contigs.fasta_name
    Description ..................................: No description is given
    Num threads for gene calling .................: 18                 

    Finding ORFs in contigs
    ===============================================
    Genes ........................................: /tmp/tmpe6hkm1vi/contigs.genes
    Amino acid sequences .........................: /tmp/tmpe6hkm1vi/contigs.amino_acid_sequences
    Log file .....................................: /tmp/tmpe6hkm1vi/00_log.txt

    CITATION
    ===============================================
    Anvi'o will use 'prodigal' by Hyatt et al (doi:10.1186/1471-2105-11-119) to
    identify open reading frames in your data. When you publish your findings,
    please do not forget to properly credit their work.

    Result .......................................: Prodigal (v2.6.3) has identified 5328 genes.

                                                                      
    CONTIGS DB CREATE REPORT
    ===============================================
    Split Length .................................: 20,000
    K-mer size ...................................: 4
    Skip gene calling? ...........................: False
    External gene calls provided? ................: False
    Ignoring internal stop codons? ...............: False
    Splitting pays attention to gene calls? ......: True
    Contigs with at least one gene call ..........: 555 of 581 (95.5%) 
    Contigs database .............................: A new database, B_theta_B_theta_WH502-contigs.fasta.db, has been created.
    Number of contigs ............................: 581
    Number of splits .............................: 839
    Total number of nucleotides ..................: 6,295,388
    Gene calling step skipped ....................: False
    Splits broke genes (non-mindful mode) ........: False
    Desired split length (what the user wanted) ..: 20,000
    Average split length (what anvi'o gave back) .: 20,968

    ✓ anvi-gen-contigs-database took 0:00:19.500103

    Working on B_theta_WH507-contigs.fasta ...

    Input FASTA file .............................: /global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/B_theta_test_2/B_theta_WH507-contigs.fasta
    Name .........................................: B_theta_B_theta_WH507-contigs.fasta_name
    Description ..................................: No description is given
    Num threads for gene calling .................: 18                 

    Finding ORFs in contigs
    ===============================================
    Genes ........................................: /tmp/tmpirdzxf3y/contigs.genes
    Amino acid sequences .........................: /tmp/tmpirdzxf3y/contigs.amino_acid_sequences
    Log file .....................................: /tmp/tmpirdzxf3y/00_log.txt

    CITATION
    ===============================================
    Anvi'o will use 'prodigal' by Hyatt et al (doi:10.1186/1471-2105-11-119) to
    identify open reading frames in your data. When you publish your findings,
    please do not forget to properly credit their work.

                                        Result .......................................: Prodigal (v2.6.3) has identified 4833 genes.

                                        [20 Aug 24 16:46:47 Processing] Enter                                                            
    CONTIGS DB CREATE REPORT
    ===============================================
    Split Length .................................: 20,000
    K-mer size ...................................: 4
    Skip gene calling? ...........................: False
    External gene calls provided? ................: False
    Ignoring internal stop codons? ...............: False
    Splitting pays attention to gene calls? ......: True
                                        [20 Aug 24 16:46:47 The South Loop] C[20 Aug 24 16:46:47 The South Loop] C[20 Aug 24 16:46:47 The South Loop] C[20 Aug 24 16:46:47 The South Loop] Cont                                        [20 Aug 24 16:46:47 The South Loop] Cont[20 Aug 24 16:46:47 The South Loop] Cont[20 Aug 24 16:46:47 The South Loop] Cont[20 Aug 24 16:46:47 The South Loop] ContContigs with at least one gene call ..........: 87 of 94 (92.6%)   
    Contigs database .............................: A new database, B_theta_B_theta_WH507-contigs.fasta.db, has been created.
    Number of contigs ............................: 94
    Number of splits .............................: 351
    Total number of nucleotides ..................: 6,157,592
    Gene calling step skipped ....................: False
    Splits broke genes (non-mindful mode) ........: False
    Desired split length (what the user wanted) ..: 20,000
    Average split length (what anvi'o gave back) .: 21,032

    ✓ anvi-gen-contigs-database took 0:00:18.159435

    Working on ref-contigs.fasta ...

    Input FASTA file .............................: /global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/B_theta_test_2/ref-contigs.fasta
    Name .........................................: B_theta_ref-contigs.fasta_name
    Description ..................................: No description is given
                                                                      

    Config Error: At least one of the deflines in your FASTA File does not comply with the 'simple
                  deflines' requirement of anvi'o. You can either use the script `anvi-script-    
                  reformat-fasta` to take care of this issue, or read this section in the tutorial
                  to understand the reason behind this requirement (anvi'o is very upset for      
                  making you do this): http://merenlab.org/2016/06/22/anvio-tutorial-v2/#take-a-  
                  look-at-your-fasta-file  

anvi-script-reformat-fasta ref-contigs.fasta --simplify-names -o ref_reformated-contigs.fasta
  Input ........................................: ref-contigs.fasta
  Output .......................................: ref_reformated-contigs.fasta

  WHAT WAS THERE
  ===============================================
  Total num contigs ............................: 2
  Total num nucleotides ........................: 6,304,193

  WHAT WAS ASKED
  ===============================================
  Simplify deflines? ...........................: Yes
  Add prefix to sequence names? ................: No
  Minimum length of contigs to keep ............: 0
  Max % gaps allowed ...........................: 100.00%
  Max num gaps allowed .........................: 1,000,000
  Exclude specific sequences? ..................: No
  Keep specific sequences? .....................: No
  Enforce sequence type? .......................: No

  WHAT HAPPENED
  ===============================================
  Contigs removed ..............................: 0 (0.00% of all)
  Nucleotides removed ..........................: 0 (0.00% of all)
  Nucleotides modified .........................: 0 (0.00000% of all)
  Deflines simplified ..........................: True


anvi-gen-contigs-database -f ref_reformated-contigs.fasta -o ref-contigs.fasta.db --num-threads 18 -n ref-contigs_name
  Input FASTA file .............................: /global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/B_theta_test_2/ref_reformated-contigs.fasta
  Name .........................................: ref-contigs_name
  Description ..................................: No description is given
  Num threads for gene calling .................: 18                               

  Finding ORFs in contigs
  ===============================================
  Genes ........................................: /tmp/tmpphgy1s00/contigs.genes
  Amino acid sequences .........................: /tmp/tmpphgy1s00/contigs.amino_acid_sequences
  Log file .....................................: /tmp/tmpphgy1s00/00_log.txt

  CITATION
  ===============================================
  Anvi'o will use 'prodigal' by Hyatt et al (doi:10.1186/1471-2105-11-119) to
  identify open reading frames in your data. When you publish your findings,
  please do not forget to properly credit their work.

                                                                                  
  WARNING
  ===============================================
  Even though you set the number of threads to 18, your FASTA file contains only 2
  sequences. To avoid any hiccups later, anvi'o will set the number of threads to
  match the number of sequences in your FASTA (who would've thought a perfect
  assembly can have a downside?).

  Result .......................................: Prodigal (v2.6.3) has identified 4945 genes.

                                                                                  
  CONTIGS DB CREATE REPORT
  ===============================================
  Split Length .................................: 20,000
  K-mer size ...................................: 4
  Skip gene calling? ...........................: False
  External gene calls provided? ................: False
  Ignoring internal stop codons? ...............: False
  Splitting pays attention to gene calls? ......: True
  Contigs with at least one gene call ..........: 2 of 2 (100.0%)                  
  Contigs database .............................: A new database, ref-contigs.fasta.db, has been created.
  Number of contigs ............................: 2
  Number of splits .............................: 314
  Total number of nucleotides ..................: 6,304,193
  Gene calling step skipped ....................: False
  Splits broke genes (non-mindful mode) ........: False
  Desired split length (what the user wanted) ..: 20,000
  Average split length (what anvi'o gave back) .: 20,036

  ✓ anvi-gen-contigs-database took 0:00:52.724719

  for g in *.db
> do
> anvi-estimate-scg-taxonomy -c $g
> done
                                                                                 
  Config Error: It seems the SCG taxonomy tables were not populated in this contigs database :/
                Luckily it is easy to fix that. Please see the program `anvi-run-scg-taxonomy`.

# so you do have to "annotate" it first

for g in *.db;
do     
anvi-run-hmms -c $g --num-threads 18;     
anvi-run-ncbi-cogs -c $g --num-threads 18 ;     
anvi-scan-trnas -c $g --num-threads 18;     
anvi-run-scg-taxonomy -c $g --num-threads 18 ; 
done

for g in *.db
do
anvi-estimate-scg-taxonomy -c $g
done
  Contigs DB ...................................: B_theta_B_theta_VPI-C11-15-contigs.fasta.db                                                      
  Metagenome mode ..............................: False
                                                                                                                                                  
  Estimated taxonomy for "B_theta_B_theta_VPI-C11-15-contigs.fasta_name"
  ===============================================
  +-----------------------------------------------+--------------+-------------------+----------------------------------------------------------------------------------------+
  |                                               |   total_scgs |   supporting_scgs | taxonomy                                                                               |
  +===============================================+==============+===================+========================================================================================+
  | B_theta_B_theta_VPI-C11-15-contigs.fasta_name |           22 |                21 | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / |
  +-----------------------------------------------+--------------+-------------------+----------------------------------------------------------------------------------------+
  Contigs DB ...................................: B_theta_B_theta_WH502-contigs.fasta.db                                                           
  Metagenome mode ..............................: False
                                                                                                                                                  
  Estimated taxonomy for "B_theta_B_theta_WH502-contigs.fasta_name"
  ===============================================
  +------------------------------------------+--------------+-------------------+----------------------------------------------------------------------------------------+
  |                                          |   total_scgs |   supporting_scgs | taxonomy                                                                               |
  +==========================================+==============+===================+========================================================================================+
  | B_theta_B_theta_WH502-contigs.fasta_name |           22 |                21 | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / |
  +------------------------------------------+--------------+-------------------+----------------------------------------------------------------------------------------+
  Contigs DB ...................................: B_theta_B_theta_WH507-contigs.fasta.db                                                                                                               
  Metagenome mode ..............................: False
                                                                                                                                                                                                      
  Estimated taxonomy for "B_theta_B_theta_WH507-contigs.fasta_name"
  ===============================================
  +------------------------------------------+--------------+-------------------+----------------------------------------------------------------------------------------+
  |                                          |   total_scgs |   supporting_scgs | taxonomy                                                                               |
  +==========================================+==============+===================+========================================================================================+
  | B_theta_B_theta_WH507-contigs.fasta_name |           22 |                21 | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / |
  +------------------------------------------+--------------+-------------------+----------------------------------------------------------------------------------------+
  Contigs DB ...................................: ref-contigs.fasta.db                                                                                                                                 
  Metagenome mode ..............................: False
                                                                                                                                                                                                      
  Estimated taxonomy for "ref-contigs_name"
  ===============================================
  +------------------+--------------+-------------------+----------------------------------------------------------------------------------------+
  |                  |   total_scgs |   supporting_scgs | taxonomy                                                                               |
  +==================+==============+===================+========================================================================================+
  | ref-contigs_name |           22 |                21 | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / |
  +------------------+--------------+-------------------+----------------------------------------------------------------------------------------+

anvi-estimate-scg-taxonomy -c ref-contigs.fasta --debug
                                                                                                                                                                                                     
  Traceback for debugging
  ================================================================================
    File "/global/home/users/empchase/.conda/envs/anvio-8-2/bin/anvi-estimate-scg-taxonomy", line 104, in <module>
      main(args)
    File "/global/home/users/empchase/.conda/envs/anvio-8-2/lib/python3.10/site-packages/anvio/terminal.py", line 915, in wrapper
      program_method(*args, **kwargs)
    File "/global/home/users/empchase/.conda/envs/anvio-8-2/bin/anvi-estimate-scg-taxonomy", line 37, in main
      t = scgtaxonomyops.SCGTaxonomyEstimatorSingle(args)
    File "/global/home/users/empchase/.conda/envs/anvio-8-2/lib/python3.10/site-packages/anvio/taxonomyops/scg.py", line 1017, in __init__
      SanityCheck.__init__(self)
    File "/global/home/users/empchase/.conda/envs/anvio-8-2/lib/python3.10/site-packages/anvio/taxonomyops/scg.py", line 135, in __init__
      self.sanity_check()
    File "/global/home/users/empchase/.conda/envs/anvio-8-2/lib/python3.10/site-packages/anvio/taxonomyops/scg.py", line 226, in sanity_check
      utils.is_contigs_db(self.contigs_db_path)
    File "/global/home/users/empchase/.conda/envs/anvio-8-2/lib/python3.10/site-packages/anvio/utils.py", line 4131, in is_contigs_db
      dbi(db_path, expecting='contigs', dont_raise=dont_raise)
    File "/global/home/users/empchase/.conda/envs/anvio-8-2/lib/python3.10/site-packages/anvio/dbinfo.py", line 93, in __new__
      if not cls.is_db(path, dont_raise=dont_raise):
    File "/global/home/users/empchase/.conda/envs/anvio-8-2/lib/python3.10/site-packages/anvio/dbinfo.py", line 170, in is_db
      raise ConfigError(f"Someone downstream doesn't like your so called database, '{path}'. They say "
  ================================================================================


    Config Error: Someone downstream doesn't like your so called database', 'ref-contigs.fasta'.
                  They say "file is not a database". Awkward :(                                


(anvio-8-2) [empchase@n0203 B_theta_test_2]$ anvi-estimate-scg-taxonomy -c ref-contigs.fasta.db --debug
Contigs DB ...................................: ref-contigs.fasta.db                                                                                                                                 
Metagenome mode ..............................: False
                                                                                                                                                                                                     
* A total of 22 single-copy core genes with taxonomic affiliations were
  successfully initialized from the contigs database 🎉 Following shows the
  frequency of these SCGs: Ribosomal_S2 (1), Ribosomal_S3_C (1), Ribosomal_S6
  (1), Ribosomal_S7 (1), Ribosomal_S8 (1), Ribosomal_S9 (1), Ribosomal_S11 (1),
  Ribosomal_S20p (1), Ribosomal_L1 (1), Ribosomal_L2 (1), Ribosomal_L3 (1),
  Ribosomal_L4 (1), Ribosomal_L6 (1), Ribosomal_L9_C (1), Ribosomal_L13 (1),
  Ribosomal_L16 (1), Ribosomal_L17 (1), Ribosomal_L20 (1), Ribosomal_L21p (1),
  Ribosomal_L22 (1), ribosomal_L24 (1), Ribosomal_L27A (1).
                                                                                                                                                                                                     
Hits for ref-contigs_name
===============================================
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| SCG            | gene   | pct id   | taxonomy                                                                                                            |
+================+========+==========+=====================================================================================================================+
| Ribosomal_L16  | 579    | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_S3_C | 580    | 99.1     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L22  | 581    | 98.1     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L2   | 583    | 98.9     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L20  | 3144   | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L4   | 585    | 99.5     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L3   | 586    | 99.0     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_S7   | 589    | 98.1     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L9_C | 18     | 99.3     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_S6   | 20     | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L1   | 596    | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_S20p | 1319   | 96.6     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / Bacteroides thetaiotaomicron |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L17  | 558    | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_S11  | 561    | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae /  /                                         |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L13  | 1778   | 98.6     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_S9   | 1779   | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / Bacteroides thetaiotaomicron |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_S2   | 1780   | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L21p | 2230   | 98.0     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / Bacteroides thetaiotaomicron |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L27A | 567    | 97.9     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L6   | 571    | 95.1     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_S8   | 572    | 97.7     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| ribosomal_L24  | 575    | 99.0     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| CONSENSUS      | --     | --       | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
                                                                                                                                                                                                     
Estimated taxonomy for "ref-contigs_name"
===============================================
+------------------+--------------+-------------------+----------------------------------------------------------------------------------------+
|                  |   total_scgs |   supporting_scgs | taxonomy                                                                               |
+==================+==============+===================+========================================================================================+
| ref-contigs_name |           22 |                21 | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / |
+------------------+--------------+-------------------+----------------------------------------------------------------------------------------+

anvi-estimate-scg-taxonomy -c B_theta_B_theta_WH502-contigs.fasta.db --debug
Contigs DB ...................................: B_theta_B_theta_WH502-contigs.fasta.db
Metagenome mode ..............................: False
                                                                                  
* A total of 22 single-copy core genes with taxonomic affiliations were
  successfully initialized from the contigs database 🎉 Following shows the
  frequency of these SCGs: Ribosomal_S2 (1), Ribosomal_S3_C (1), Ribosomal_S6
  (1), Ribosomal_S7 (1), Ribosomal_S8 (1), Ribosomal_S9 (1), Ribosomal_S11 (1),
  Ribosomal_S20p (1), Ribosomal_L1 (1), Ribosomal_L2 (1), Ribosomal_L3 (1),
  Ribosomal_L4 (1), Ribosomal_L6 (1), Ribosomal_L9_C (1), Ribosomal_L13 (1),
  Ribosomal_L16 (1), Ribosomal_L17 (1), Ribosomal_L20 (1), Ribosomal_L21p (1),
  Ribosomal_L22 (1), ribosomal_L24 (1), Ribosomal_L27A (1).
                                                                                  
Hits for B_theta_B_theta_WH502-contigs.fasta_name
===============================================
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| SCG            | gene   | pct id   | taxonomy                                                                                                  |
+================+========+==========+===========================================================================================================+
| Ribosomal_S9   | 576    | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / Bacteroides faecis |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L13  | 577    | 98.6     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L27A | 3521   | 97.2     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_S11  | 3527   | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae /  /                               |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L17  | 3530   | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L9_C | 1549   | 99.3     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / Bacteroides faecis |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_S6   | 1551   | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_S20p | 4573   | 96.6     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L1   | 3492   | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L20  | 2726   | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_S7   | 3499   | 98.1     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L3   | 3502   | 99.0     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L4   | 3503   | 98.4     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L2   | 3505   | 98.5     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L21p | 115    | 98.0     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / Bacteroides faecis |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L22  | 3507   | 98.1     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_S3_C | 3508   | 99.1     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L16  | 3509   | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| ribosomal_L24  | 3513   | 98.1     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_S8   | 3516   | 97.7     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L6   | 3517   | 95.1     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_S2   | 575    | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| CONSENSUS      | --     | --       | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
                                                                                  
Estimated taxonomy for "B_theta_B_theta_WH502-contigs.fasta_name"
===============================================
+------------------------------------------+--------------+-------------------+----------------------------------------------------------------------------------------+
|                                          |   total_scgs |   supporting_scgs | taxonomy                                                                               |
+==========================================+==============+===================+========================================================================================+
| B_theta_B_theta_WH502-contigs.fasta_name |           22 |                21 | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / |
+------------------------------------------+--------------+-------------------+----------------------------------------------------------------------------------------+
(anvio-8-2) [empchase@n0203 B_theta_test_2]$ anvi-estimate-scg-taxonomy -c B_theta_B_theta_WH507-contigs.fasta.db --debug
Contigs DB ...................................: B_theta_B_theta_WH507-contigs.fasta.db                                                                                                               
Metagenome mode ..............................: False
                                                                                                                                                                                                     
* A total of 22 single-copy core genes with taxonomic affiliations were
  successfully initialized from the contigs database 🎉 Following shows the
  frequency of these SCGs: Ribosomal_S2 (1), Ribosomal_S3_C (1), Ribosomal_S6
  (1), Ribosomal_S7 (1), Ribosomal_S8 (1), Ribosomal_S9 (1), Ribosomal_S11 (1),
  Ribosomal_S20p (1), Ribosomal_L1 (1), Ribosomal_L2 (1), Ribosomal_L3 (1),
  Ribosomal_L4 (1), Ribosomal_L6 (1), Ribosomal_L9_C (1), Ribosomal_L13 (1),
  Ribosomal_L16 (1), Ribosomal_L17 (1), Ribosomal_L20 (1), Ribosomal_L21p (1),
  Ribosomal_L22 (1), ribosomal_L24 (1), Ribosomal_L27A (1).
                                                                                                                                                                                                     
Hits for B_theta_B_theta_WH507-contigs.fasta_name
===============================================
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| SCG            | gene   | pct id   | taxonomy                                                                                                  |
+================+========+==========+===========================================================================================================+
| Ribosomal_S11  | 3330   | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae /  /                               |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L27A | 3336   | 97.2     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L20  | 2827   | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L6   | 3340   | 95.1     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_S8   | 3341   | 97.7     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| ribosomal_L24  | 3344   | 98.1     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L16  | 3348   | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_S3_C | 3349   | 99.1     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L22  | 3350   | 98.1     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L2   | 3352   | 98.5     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L4   | 3354   | 98.4     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L3   | 3355   | 99.0     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_S7   | 3358   | 98.1     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L9_C | 2208   | 99.3     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / Bacteroides faecis |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_S6   | 2210   | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L1   | 3365   | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L21p | 3944   | 98.0     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / Bacteroides faecis |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L13  | 42     | 98.6     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_S9   | 43     | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / Bacteroides faecis |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_S2   | 44     | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_S20p | 4541   | 96.6     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| Ribosomal_L17  | 3327   | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
| CONSENSUS      | --     | --       | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                    |
+----------------+--------+----------+-----------------------------------------------------------------------------------------------------------+
                                                                                                                                                                                                     
Estimated taxonomy for "B_theta_B_theta_WH507-contigs.fasta_name"
===============================================
+------------------------------------------+--------------+-------------------+----------------------------------------------------------------------------------------+
|                                          |   total_scgs |   supporting_scgs | taxonomy                                                                               |
+==========================================+==============+===================+========================================================================================+
| B_theta_B_theta_WH507-contigs.fasta_name |           22 |                21 | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / |
+------------------------------------------+--------------+-------------------+----------------------------------------------------------------------------------------+
(anvio-8-2) [empchase@n0203 B_theta_test_2]$ anvi-estimate-scg-taxonomy -c B_theta_VPI-C11-15-contigs.fasta --debug
B_theta_VPI-C11-15-contigs.fasta
(anvio-8-2) [empchase@n0203 B_theta_test_2]$ anvi-estimate-scg-taxonomy -c B_theta_B_theta_VPI-C11-15-contigs.fasta.db --debug
Contigs DB ...................................: B_theta_B_theta_VPI-C11-15-contigs.fasta.db                                                                                                          
Metagenome mode ..............................: False
                                                                                                                                                                                                     
* A total of 22 single-copy core genes with taxonomic affiliations were
  successfully initialized from the contigs database 🎉 Following shows the
  frequency of these SCGs: Ribosomal_S2 (1), Ribosomal_S3_C (1), Ribosomal_S6
  (1), Ribosomal_S7 (1), Ribosomal_S8 (1), Ribosomal_S9 (1), Ribosomal_S11 (1),
  Ribosomal_S20p (1), Ribosomal_L1 (1), Ribosomal_L2 (1), Ribosomal_L3 (1),
  Ribosomal_L4 (1), Ribosomal_L6 (1), Ribosomal_L9_C (1), Ribosomal_L13 (1),
  Ribosomal_L16 (1), Ribosomal_L17 (1), Ribosomal_L20 (1), Ribosomal_L21p (1),
  Ribosomal_L22 (1), ribosomal_L24 (1), Ribosomal_L27A (1).
                                                                                                                                                                                                     
Hits for B_theta_B_theta_VPI-C11-15-contigs.fasta_name
===============================================
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| SCG            | gene   | pct id   | taxonomy                                                                                                            |
+================+========+==========+=====================================================================================================================+
| Ribosomal_S6   | 1544   | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_S2   | 3594   | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L9_C | 1546   | 99.3     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_S9   | 3595   | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / Bacteroides faecis           |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L13  | 3596   | 98.6     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L1   | 5518   | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_S7   | 5525   | 98.1     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L3   | 5528   | 99.0     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L4   | 5529   | 99.5     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L2   | 5531   | 98.9     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L22  | 5533   | 98.1     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_S3_C | 5534   | 99.1     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L16  | 5535   | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| ribosomal_L24  | 5539   | 99.0     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_S8   | 5542   | 97.7     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L6   | 5543   | 95.1     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L20  | 3687   | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L21p | 5290   | 98.0     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / Bacteroides thetaiotaomicron |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L27A | 5547   | 97.9     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_S20p | 3184   | 96.6     | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / Bacteroides thetaiotaomicron |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_S11  | 5553   | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae /  /                                         |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| Ribosomal_L17  | 5556   | 100.0    | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
| CONSENSUS      | --     | --       | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides /                              |
+----------------+--------+----------+---------------------------------------------------------------------------------------------------------------------+
                                                                                                                                                                                                     
Estimated taxonomy for "B_theta_B_theta_VPI-C11-15-contigs.fasta_name"
===============================================
+-----------------------------------------------+--------------+-------------------+----------------------------------------------------------------------------------------+
|                                               |   total_scgs |   supporting_scgs | taxonomy                                                                               |
+===============================================+==============+===================+========================================================================================+
| B_theta_B_theta_VPI-C11-15-contigs.fasta_name |           22 |                21 | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / |
+-----------------------------------------------+--------------+-------------------+----------------------------------------------------------------------------------------+

# August 21, 2024
# How can we detect contamination?

#checkm on 3 bthetas and reference genome
cd /global/scratch/projects/fc_wolflab/bbioinfo/empchase/checkm_test
# rewrote checkm_test.sh to look at all fasta files in /global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/B_theta_test_2
sbatch checkm_test.sh 
  Submitted batch job 20569471



#does anvio pangenomics agree?
cd /global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/B_theta_test_2
# waiting on savio interactive




# checkm on synthetic contaminants
# step1 create synthetic contaminant file: 2 b theta strains
# Going to use A8 and A4 B theta's bc they have 88% ANI
cd /global/scratch/projects/fc_wolflab/bbioinfo/empchase/synthetic_contaminant

#organize
[empchase@n0000 synthetic_contaminant]$ mkdir D4theta90_F1vulgatus10
[empchase@n0000 synthetic_contaminant]$ mv *.gz D4theta90_F1vulgatus10/
[empchase@n0000 synthetic_contaminant]$ mv spades_test/ D4theta90_F1vulgatus10/
[empchase@n0000 synthetic_contaminant]$ mv test/ D4theta90_F1vulgatus10/
[empchase@n0000 synthetic_contaminant]$ mv trimmomatic_test/ D4theta90_F1vulgatus10/
[empchase@n0000 synthetic_contaminant]$ mv reads_analysis.ipynb D4theta90_F1vulgatus10/

[empchase@n0000 synthetic_contaminant]$ mkdir A8thetaVPI-90_A4thetaWH502-10
#edit contaminate_v2.sh and contaminate_v2.py to reflect file paths
conda activate cmd_biopython
sbatch contaminate_v2.sh 
  Submitted batch job 20569568

# August 22, 2024

# contaminant
# moved all subsetA4/A8 files to /global/scratch/projects/fc_wolflab/bbioinfo/empchase/synthetic_contaminant/A8thetaVPI-90_A4thetaWH502-10

#trimmomatic and spades
cat /global/scratch/projects/fc_wolflab/bbioinfo/sarrah/batches/adapter_files/A8_adapters.fasta /global/scratch/projects/fc_wolflab/bbioinfo/sarrah/batches/adapter_files/A4_adapters.fasta > A8+A4_adapters.fasta
conda activate python3-8-8
sbatch trim_spades.sh 
  Submitted batch job 20597666

# anvio
conda deactivate
conda activate anvio-8-2

anvi-script-gen-genomes-file --input-dir . -o external-genomes.txt

anvi-gen-genomes-storage -e external-genomes.txt -o $B_theta-GENOMES.db
  usage: anvi-gen-genomes-storage [-h] [-e FILE_PATH] [-i FILE_PATH] [--gene-caller GENE-CALLER] -o GENOMES_STORAGE
  anvi-gen-genomes-storage: error: argument -o/--output-file: expected one argument

anvi-gen-genomes-storage -e external-genomes.txt -o B_theta-GENOMES.db
  Config Error: Well, the genome name "B_theta_B_theta_VPI-C11-15-contigs.fasta_name" contains
                characters that anvi'o does not like :/ Please limit the characters to ASCII  
                letters, digits, and the underscore ('_') character.   '


(anvio-8-2) [empchase@n0206 B_theta_test_2]$ mkdir oops
(anvio-8-2) [empchase@n0206 B_theta_test_2]$ mv *.db oops
(anvio-8-2) [empchase@n0206 B_theta_test_2]$ ls
B_theta_VPI-C11-15-contigs.fasta  B_theta_WH502-contigs.fasta  B_theta_WH507-contigs.fasta  external-genomes.txt  oops  ref-contigs.fasta  ref_reformated-contigs.fasta
(anvio-8-2) [empchase@n0206 B_theta_test_2]$ mv B_theta_VPI-C11-15-contigs.fasta B_theta_VPIC1115_contigs.fasta
(anvio-8-2) [empchase@n0206 B_theta_test_2]$ mv B_theta_WH502-contigs.fasta B_theta_WH502_contigs.fasta
(anvio-8-2) [empchase@n0206 B_theta_test_2]$ mv B_theta_WH507-contigs.fasta B_theta_WH507_contigs.fasta
(anvio-8-2) [empchase@n0206 B_theta_test_2]$ ls
B_theta_VPIC1115_contigs.fasta  B_theta_WH502_contigs.fasta  B_theta_WH507_contigs.fasta  external-genomes.txt  oops  ref-contigs.fasta  ref_reformated-contigs.fasta
(anvio-8-2) [empchase@n0206 B_theta_test_2]$ mv external-genomes.txt oops
(anvio-8-2) [empchase@n0206 B_theta_test_2]$ mv ref_reformated-contigs.fasta ref_reformated_contigs.fasta
(anvio-8-2) [empchase@n0206 B_theta_test_2]$ for g in *.fasta
> do
> anvi-gen-contigs-database -f $g -o ${g}.db --num-threads 20 -n ${g}_name
> done
  
for g in *.db;
do     
anvi-run-hmms -c $g --num-threads 18;     
anvi-run-ncbi-cogs -c $g --num-threads 18 ;     
anvi-scan-trnas -c $g --num-threads 18;     
anvi-run-scg-taxonomy -c $g --num-threads 18 ; 
done

anvi-script-gen-genomes-file --input-dir . -o external-genomes.txt
anvi-gen-genomes-storage -e external-genomes.txt -o B_theta-GENOMES.db
  Config Error: Well, the genome name "B_theta_VPIC1115_contigs.fasta_name" contains characters
                that anvi'o does not like :/ Please limit the characters to ASCII letters,     
                digits, and the underscore ('_') character. '
                
# anvi-pan-genome -g B_theta-GENOMES.db --project-name B_theta --num-threads 18 # not run yet

# August 26, 2024
# checkm btheta
sbatch /global/scratch/projects/fc_wolflab/bbioinfo/empchase/checkm_test/checkm_test.sh
  Submitted batch job 20623008

# anvio fix genome name
cd /global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/B_theta_test_2
# just gonna redo again ugh
#reorganize
[empchase@n0087 B_theta_test_2]$ mv *.db oops
[empchase@n0087 B_theta_test_2]$ ls
B_theta_VPIC1115_contigs.fasta  B_theta_WH502_contigs.fasta  B_theta_WH507_contigs.fasta  external-genomes.txt  oops  ref-contigs.fasta  ref_reformated_contigs.fasta
[empchase@n0087 B_theta_test_2]$ mv external-genomes.txt external-genomes2.txt
[empchase@n0087 B_theta_test_2]$ mv external-genomes2.txt oops/
[empchase@n0087 B_theta_test_2]$ ls
B_theta_VPIC1115_contigs.fasta  B_theta_WH502_contigs.fasta  B_theta_WH507_contigs.fasta  oops  ref-contigs.fasta  ref_reformated_contigs.fasta
[empchase@n0087 B_theta_test_2]$ mkdir old
[empchase@n0087 B_theta_test_2]$ mv ref-contigs.fasta old

#redo
ls *.fasta | cut -d "." -f 1 > strain_names.txt

# for g in `cat strain_names.txt`
# do
# head -n 2 ${g}.fasta 
# done

for g in `cat strain_names.txt`
do
anvi-gen-contigs-database -f ${g}.fasta -o ${g}.db --num-threads 20 -n ${g}_name
done
  
for g in *.db;
do     
anvi-run-hmms -c $g --num-threads 20;     
anvi-run-ncbi-cogs -c $g --num-threads 20 ;     
anvi-scan-trnas -c $g --num-threads 20;     
anvi-run-scg-taxonomy -c $g --num-threads 20 ; 
done

anvi-script-gen-genomes-file --input-dir . -o external-genomes.txt
anvi-gen-genomes-storage -e external-genomes.txt -o B_theta-GENOMES.db
anvi-pan-genome -g B_theta-GENOMES.db --project-name B_theta --num-threads 20 # not run yet

anvi-compute-genome-similarity --external-genomes external-genomes.txt --program pyANI --output-dir ANI --num-threads 20 --pan-db B_theta/B_theta-PAN.db

for g in *_contigs.db
> do
> anvi-estimate-scg-taxonomy -c $g
> done
Contigs DB ...................................: B_theta_VPIC1115_contigs.db            
Metagenome mode ..............................: False
                                                                                       
Estimated taxonomy for "B_theta_VPIC1115_contigs_name"
===============================================
+-------------------------------+--------------+-------------------+----------------------------------------------------------------------------------------+
|                               |   total_scgs |   supporting_scgs | taxonomy                                                                               |
+===============================+==============+===================+========================================================================================+
| B_theta_VPIC1115_contigs_name |           22 |                21 | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / |
+-------------------------------+--------------+-------------------+----------------------------------------------------------------------------------------+
Contigs DB ...................................: B_theta_WH502_contigs.db               
Metagenome mode ..............................: False
                                                                                       
Estimated taxonomy for "B_theta_WH502_contigs_name"
===============================================
+----------------------------+--------------+-------------------+----------------------------------------------------------------------------------------+
|                            |   total_scgs |   supporting_scgs | taxonomy                                                                               |
+============================+==============+===================+========================================================================================+
| B_theta_WH502_contigs_name |           22 |                21 | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / |
+----------------------------+--------------+-------------------+----------------------------------------------------------------------------------------+
Contigs DB ...................................: B_theta_WH507_contigs.db               
Metagenome mode ..............................: False
                                                                                                                                                                             
Estimated taxonomy for "B_theta_WH507_contigs_name"
===============================================
+----------------------------+--------------+-------------------+----------------------------------------------------------------------------------------+
|                            |   total_scgs |   supporting_scgs | taxonomy                                                                               |
+============================+==============+===================+========================================================================================+
| B_theta_WH507_contigs_name |           22 |                21 | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / |
+----------------------------+--------------+-------------------+----------------------------------------------------------------------------------------+
Contigs DB ...................................: ref_reformated_contigs.db                                                                                                    
Metagenome mode ..............................: False
                                                                                                                                                                             
Estimated taxonomy for "ref_reformated_contigs_name"
===============================================
+-----------------------------+--------------+-------------------+----------------------------------------------------------------------------------------+
|                             |   total_scgs |   supporting_scgs | taxonomy                                                                               |
+=============================+==============+===================+========================================================================================+
| ref_reformated_contigs_name |           22 |                21 | Bacteria / Bacteroidota / Bacteroidia / Bacteroidales / Bacteroidaceae / Bacteroides / |
+-----------------------------+--------------+-------------------+----------------------------------------------------------------------------------------+

# August 27, 2024
cd /global/scratch/projects/fc_wolflab/bbioinfo/empchase/anvio_test/B_theta_test_2
conda activate anvio-8-2

for g in *_contigs.db
do
anvi-estimate-scg-taxonomy -c $g --debug
done
# see notebook/slack for results

#checkm on all genomes
cd /global/scratch/projects/fc_wolflab/bbioinfo/empchase
mkdir checkm_all
cd checkm_all

#copy /global/scratch/projects/fc_wolflab/bbioinfo/empchase/checkm_test/checkm_test.sh 
touch checkm.sh #mkdir results
# added options to multithread
sbatch checkm.sh
  Submitted batch job 20668489
# ran successfully, analyzing via jupyter notebook
# pasted table from slurm.out into outputtable.txt
# paste table from Bacteroides StrainsID.xlsx from Vienvilay/Alyssa into 16Sconfirmed_20240827.tsv

# August 29, 2024
# checkm and contamination detection

# run checkm on /global/scratch/projects/fc_wolflab/bbioinfo/empchase/synthetic_contaminant/A8thetaVPI-90_A4thetaWH502-10
# edit /global/scratch/projects/fc_wolflab/bbioinfo/empchase/checkm_test/checkm_contaminant.sh
# input_fastadir="/global/scratch/projects/fc_wolflab/bbioinfo/empchase/synthetic_contaminant/A8thetaVPI-90_A4thetaWH502-10"
# output_dir="/global/scratch/projects/fc_wolflab/bbioinfo/empchase/checkm_test/contamBtheta_A8vA4_results_rocky8"
cd /global/scratch/projects/fc_wolflab/bbioinfo/empchase/checkm_test
sbatch checkm_contaminant.sh 
  Submitted batch job 20706673 # error bc adaptors fasta file was in the directory
mv /global/scratch/projects/fc_wolflab/bbioinfo/empchase/synthetic_contaminant/A8thetaVPI-90_A4thetaWH502-10/A8+A4_adapters.fasta /global/scratch/projects/fc_wolflab/bbioinfo/empchase/synthetic_contaminant/A8thetaVPI-90_A4thetaWH502-10/trimmomatic_test/
sbatch checkm_contaminant.sh
  Submitted batch job 20707404

# August 30, 2024
# retry generating synthetic contam 
# create synthetic contaminant from 2 verified strains w/in a species and do the same workflow (spades, checkm)
  
# A02 B. vulgatus	WH019
# A08 B. vulgatus	WH515

# September 18, 2024
# I want to know how to interpret the kraken data sarrah generated found in /global/scratch/projects/fc_wolflab/bbioinfo/sarrah/batches/trimmed_all/
mkdir PCR_contamination_test
touch PCR_results.txt #tsv; similar to checkm_all/16sdata.txt but I changed the well IDs from eg A01 to A1
# get filepaths from /global/scratch/projects/fc_wolflab/bbioinfo/sarrah/batches/Corrected_Final_Config.csv, combine in PCR_filepath_generate.ipynb
# PCRtest_tilepaths.txt

mkdir readsforbvbrc
cd readsforbvbrc
for f in `cat ../PCRtest_filepaths.txt`; do cp $f .; done
# from local: scp -r empchase@dtn.brc.berkeley.edu:/global/scratch/projects/fc_wolflab/bbioinfo/empchase/PCR_contamination_test/readsforbvbrc ~/Desktop/WolfLab/Bacteroides_Bioinformatics/rawreads/
# from local: scp empchase@dtn.brc.berkeley.edu:/global/scratch/projects/fc_wolflab/bbioinfo/empchase/PCR_contamination_test/PCRtest_filepaths.txt ~/Desktop/WolfLab/Bacteroides_Bioinformatics/rawreads/readsforbvbrc


# September 19, 2024
cd /global/scratch/projects/fc_wolflab/bbioinfo/empchase/synthetic_contaminant
mkdir vulgatus_ncbi
cd vulgatus_ncbi
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/012/825/GCF_000012825.1_ASM1282v1/GCF_000012825.1_ASM1282v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/273/295/GCF_000273295.1_Bact_vulg_CL09T03C04_V1/GCF_000273295.1_Bact_vulg_CL09T03C04_V1_genomic.fna.gz