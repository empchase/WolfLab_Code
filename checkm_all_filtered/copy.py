import os
import shutil
import csv
import glob

# Path to your input CSV file
input_file = '/global/scratch/projects/fc_wolflab/bbioinfo/empchase/PCR_contamination_test/genomes_singlekrakenassignment_humanfil.csv'

# Destination directory
destination_dir = 'tmp/'

# Read the CSV file
with open(input_file, 'r') as file:
    reader = csv.DictReader(file)

    for row in reader:
        # Extract necessary information
        plate_number = row['Plate_number']
        library_name = row['LibraryName']
        genome_path = row['genome_path']
        # print(genome_path)

        # Use glob to find the actual directory that matches the wildcard pattern
        assembly_dir_pattern = os.path.join(genome_path, '*_assembly')
        assembly_dirs = glob.glob(assembly_dir_pattern)

        # Ensure that we found a matching directory
        if assembly_dirs:
            # print(len(assembly_dirs))
            # There may be multiple matches, so we'll just take the first one
            assembly_dir = assembly_dirs[0]

            # Define the source file (contigs.fasta) within the matched directory
            source_file = os.path.join(assembly_dir, 'contigs.fasta')

            # Define the destination file name with the required format
            destination_file = os.path.join(destination_dir, f'{plate_number}_{library_name}_contigs.fasta')

            # Check if the source file exists and copy it
            if os.path.exists(source_file):
                try:
                    shutil.copy(source_file, destination_file)
                    print(f'Successfully copied: {source_file} to {destination_file}')
                except Exception as e:
                    print(f'Error copying {source_file} to {destination_file}: {e}')
            else:
                print(f'contigs.fasta not found in: {assembly_dir}')
        else:
            print(f'No assembly directory found for pattern: {assembly_dir_pattern}')
