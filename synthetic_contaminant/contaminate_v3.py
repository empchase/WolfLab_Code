import argparse
from Bio import SeqIO
import random
import gzip

def process_sample(fastq1, fastq2, output1, output2, proportion):
    read1_ids = set()
    total_reads = 0

    # Process read1 file
    with gzip.open(fastq1, "rt") as handle1, gzip.open(output1, "wt") as out_handle1:
        for record in SeqIO.parse(handle1, "fastq"):
            total_reads += 1
            if random.random() < proportion:
                read1_ids.add(record.id)
                SeqIO.write(record, out_handle1, "fastq")

    # Process read2 file
    with gzip.open(fastq2, "rt") as handle2, gzip.open(output2, "wt") as out_handle2:
        for record in SeqIO.parse(handle2, "fastq"):
            if record.id in read1_ids:
                SeqIO.write(record, out_handle2, "fastq")

    return total_reads

def process_sample_with_exact_count(fastq1, fastq2, output1, output2, num_reads):
    read1_ids = set()
    total_reads = 0

    # Process read1 file
    with gzip.open(fastq1, "rt") as handle1, gzip.open(output1, "wt") as out_handle1:
        for record in SeqIO.parse(handle1, "fastq"):
            total_reads += 1
            if total_reads <= num_reads:
                read1_ids.add(record.id)
                SeqIO.write(record, out_handle1, "fastq")
            if total_reads >= num_reads:
                break

    # Process read2 file
    with gzip.open(fastq2, "rt") as handle2, gzip.open(output2, "wt") as out_handle2:
        for record in SeqIO.parse(handle2, "fastq"):
            if record.id in read1_ids:
                SeqIO.write(record, out_handle2, "fastq")
            if len(read1_ids) >= num_reads:
                break

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Process paired-end FASTQ files to create subsets with specified proportions")

    # Arguments for Sample A and Sample B inputs and outputs
    parser.add_argument("fastq1_sampleA", type=str, help="Path to Sample A Read 1 FASTQ file (gzipped)")
    parser.add_argument("fastq2_sampleA", type=str, help="Path to Sample A Read 2 FASTQ file (gzipped)")
    parser.add_argument("output_subsetA_read1", type=str, help="Output path for Sample A subset Read 1 FASTQ file (gzipped)")
    parser.add_argument("output_subsetA_read2", type=str, help="Output path for Sample A subset Read 2 FASTQ file (gzipped)")
    parser.add_argument("--proportion_A", type=float, default=0.9, help="Proportion of reads to keep from Sample A (default: 0.9)")

    parser.add_argument("fastq1_sampleB", type=str, help="Path to Sample B Read 1 FASTQ file (gzipped)")
    parser.add_argument("fastq2_sampleB", type=str, help="Path to Sample B Read 2 FASTQ file (gzipped)")
    parser.add_argument("output_subsetB_read1", type=str, help="Output path for Sample B subset Read 1 FASTQ file (gzipped)")
    parser.add_argument("output_subsetB_read2", type=str, help="Output path for Sample B subset Read 2 FASTQ file (gzipped)")

    # Parse the arguments
    args = parser.parse_args()

    # Process Sample A
    total_reads_A = process_sample(args.fastq1_sampleA, args.fastq2_sampleA, args.output_subsetA_read1, args.output_subsetA_read2, args.proportion_A)

    # Calculate the number of reads to sample from Sample B
    num_reads_B = int((1 - args.proportion_A) * total_reads_A)

    # Process Sample B with the exact number of reads
    process_sample_with_exact_count(args.fastq1_sampleB, args.fastq2_sampleB, args.output_subsetB_read1, args.output_subsetB_read2, num_reads_B)
