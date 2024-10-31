import argparse
from Bio import SeqIO
import gzip

def generate_pseudo_reads(input_fasta, output_fasta, read_length=100, step_size=50):
    # Open .gz files if needed
    if input_fasta.endswith(".gz"):
        with gzip.open(input_fasta, "rt") as infile, open(output_fasta, "w") as out_f:
            process_fasta(infile, out_f, read_length, step_size)
    else:
        with open(input_fasta, "r") as infile, open(output_fasta, "w") as out_f:
            process_fasta(infile, out_f, read_length, step_size)

def process_fasta(infile, out_f, read_length, step_size):
    for record in SeqIO.parse(infile, "fasta"):
        genome_seq = str(record.seq)  # Convert the sequence to a string
        genome_length = len(genome_seq)
        
        # Generate overlapping pseudo-reads
        read_count = 1
        for i in range(0, genome_length - read_length + 1, step_size):
            pseudo_read = genome_seq[i:i + read_length]
            
            # Write each pseudo-read to the output file
            out_f.write(f">pseudo_read_{read_count}_{i+1}-{i+read_length}\n")
            out_f.write(f"{pseudo_read}\n")
            read_count += 1

        print(f"Generated {read_count - 1} pseudo-reads from the genome.")

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Generate overlapping pseudo-reads from a genome FASTA file")
    
    # Required arguments
    parser.add_argument("input_fasta", type=str, help="Input genome FASTA file (can be .gz)")
    parser.add_argument("output_fasta", type=str, help="Output FASTA file for pseudo-reads")
    
    # Optional arguments
    parser.add_argument("--read_length", type=int, default=100, help="Length of each pseudo-read (default: 100)")
    parser.add_argument("--step_size", type=int, default=50, help="Step size between the start of each pseudo-read (default: 50)")

    # Parse arguments
    args = parser.parse_args()

    # Call the function with parsed arguments
    generate_pseudo_reads(args.input_fasta, args.output_fasta, args.read_length, args.step_size)
