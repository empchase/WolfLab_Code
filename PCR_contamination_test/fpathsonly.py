import csv

# Input and output file paths
input_file = 'genomes_singlekrakenassignment_humanfil.csv'  # Replace with your input file path
output_file = 'genomes_singlekrakenassignment_humanfil_fpathsonly.csv'  # Replace with your desired output file path

# Open the input file
with open(input_file, mode='r') as infile, open(output_file, mode='w') as outfile:
    # Read the CSV file
    reader = csv.DictReader(infile)
    
    # Iterate over each row
    for row in reader:
        # Get the genome_path and append "/genome_assembly"
        modified_path = row['genome_path'] + '/genome_assembly'
        
        # Write the modified path to the output file
        outfile.write(modified_path + '\n')

print("Modified genome paths have been saved to", output_file)
