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

# Define sample file paths
fastq1_sampleA_read1="/global/scratch/projects/fc_wolflab/bbioinfo/230914_BacteroidesPatnodeLab_GTAC_VPL/BacteroidesPatnodePlate1_121323/LIB015015_22GMJHLT3_S95_L006_R1_001.fastq.gz"
fastq2_sampleA_read2="/global/scratch/projects/fc_wolflab/bbioinfo/230914_BacteroidesPatnodeLab_GTAC_VPL/BacteroidesPatnodePlate1_121323/LIB015015_22GMJHLT3_S95_L006_R2_001.fastq.gz"
fastq1_sampleB_read1="/global/scratch/projects/fc_wolflab/bbioinfo/230914_BacteroidesPatnodeLab_GTAC_VPL/BacteroidesPatnodePlate1_121323/LIB015018_22GMJHLT3_S92_L006_R1_001.fastq.gz"
fastq2_sampleB_read2="/global/scratch/projects/fc_wolflab/bbioinfo/230914_BacteroidesPatnodeLab_GTAC_VPL/BacteroidesPatnodePlate1_121323/LIB015018_22GMJHLT3_S92_L006_R2_001.fastq.gz"
output_subsetA_read1="subsetA02_read1.fastq.gz"
output_subsetA_read2="subsetA02_read2.fastq.gz"
output_subsetB_read1="subsetA05_read1.fastq.gz"
output_subsetB_read2="subsetA05_read2.fastq.gz"

# Process Sample A: 90%
total_reads_A = process_sample(fastq1_sampleA_read1, fastq2_sampleA_read2, output_subsetA_read1, output_subsetA_read2, 0.9)

# Calculate the number of reads to sample from Sample B
num_reads_B = int(0.1 * total_reads_A)

# Process Sample B: equivalent to 10% of Sample A reads
# Note: proportion for Sample B is calculated based on the number of reads needed
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

process_sample_with_exact_count(fastq1_sampleB_read1, fastq2_sampleB_read2, output_subsetB_read1, output_subsetB_read2, num_reads_B)
