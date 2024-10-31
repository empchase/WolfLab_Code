import argparse
from Bio import SeqIO
import random
import gzip

def process_sample(fasta1, output1, proportion):
    read1_ids = set()
    total_reads = 0

    # Process read1 file
    with open(fasta1, "r") as handle1, open(output1, "w") as out_handle1:
        for record in SeqIO.parse(handle1, "fasta"):
            total_reads += 1
            if random.random() < proportion:
                read1_ids.add(record.id)
                SeqIO.write(record, out_handle1, "fasta")

    return total_reads

def process_sample_with_exact_count(fasta1, output1, num_reads):
    read1_ids = set()
    total_reads = 0

    # Process read1 file
    with open(fasta1, "r") as handle1, open(output1, "w") as out_handle1:
        for record in SeqIO.parse(handle1, "fasta"):
            total_reads += 1
            if total_reads <= num_reads:
                read1_ids.add(record.id)
                SeqIO.write(record, out_handle1, "fasta")
            if total_reads >= num_reads:
                break


if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Process FASTA files to create subsets with specified proportions")

    # Arguments for Sample A and Sample B inputs and outputs
    parser.add_argument("fasta1_sampleA", type=str, help="Path to Sample A Read 1 FASTA file (gzipped)")
    parser.add_argument("output_subsetA_read1", type=str, help="Output path for Sample A subset Read 1 FASTA file (gzipped)")
    parser.add_argument("--proportion_A", type=float, default=0.9, help="Proportion of reads to keep from Sample A (default: 0.9)")

    parser.add_argument("fasta1_sampleB", type=str, help="Path to Sample B Read 1 FASTA file (gzipped)")
    parser.add_argument("output_subsetB_read1", type=str, help="Output path for Sample B subset Read 1 FASTA file (gzipped)")

    # Parse the arguments
    args = parser.parse_args()

    # Process Sample A
    total_reads_A = process_sample(args.fasta1_sampleA, args.output_subsetA_read1, args.proportion_A)

    # Calculate the number of reads to sample from Sample B
    num_reads_B = int((1 - args.proportion_A) * total_reads_A)

    # Process Sample B with the exact number of reads
    process_sample_with_exact_count(args.fasta1_sampleB, args.output_subsetB_read1, num_reads_B)
