from Bio import SeqIO
import random
import gzip

def process_sample(fastq1, fastq2, output1, output2, proportion):
    read1_ids = set()
    total_reads = 0
    
    # Step 1-2: Process read1 file
    with gzip.open(fastq1, "rt") as handle1, gzip.open(output1, "wt") as out_handle1:
        for record in SeqIO.parse(handle1, "fastq"):
            total_reads += 1
            if random.random() < proportion:
                read1_ids.add(record.id)
                SeqIO.write(record, out_handle1, "fastq")
    
    # Step 3-4: Process read2 file
    with gzip.open(fastq2, "rt") as handle2, gzip.open(output2, "wt") as out_handle2:
        for record in SeqIO.parse(handle2, "fastq"):
            if record.id in read1_ids:
                SeqIO.write(record, out_handle2, "fastq")
    
    return total_reads

# Define sample file paths
fastq1_sampleA_read1 = "/global/scratch/projects/fc_wolflab/bbioinfo/230914_BacteroidesPatnodeLab_GTAC_VPL/BacteroidesPatnodePlate1_121323/LIB015053_22GMJHLT3_S57_L006_R1_001.fastq.gz"
fastq2_sampleA_read2 = "/global/scratch/projects/fc_wolflab/bbioinfo/230914_BacteroidesPatnodeLab_GTAC_VPL/BacteroidesPatnodePlate1_121323/LIB015053_22GMJHLT3_S57_L006_R2_001.fastq.gz"
# fastq1_sampleB_read1 = "/global/scratch/projects/fc_wolflab/bbioinfo/230914_BacteroidesPatnodeLab_GTAC_VPL/BacteroidesPatnodePlate1_121323/LIB015074_22GMJHLT3_S36_L006_R1_001.fastq.gz"
# fastq2_sampleB_read2 = "/global/scratch/projects/fc_wolflab/bbioinfo/230914_BacteroidesPatnodeLab_GTAC_VPL/BacteroidesPatnodePlate1_121323/LIB015074_22GMJHLT3_S36_L006_R2_001.fastq.gz"
output_subsetA_read1 = "subsetD4_read1.fastq.gz"
output_subsetA_read2 = "subsetD4_read2.fastq.gz"
# output_subsetB_read1 = "subsetF1_read1.fastq.gz"
# output_subsetB_read2 = "subsetF1_read2.fastq.gz"

# Process Sample A: 90%
total_reads_A = process_sample(fastq1_sampleA_read1, fastq2_sampleA_read2, output_subsetA_read1, output_subsetA_read2, 0.9)

# Calculate 10% of total reads from Sample A for Sample B
proportion_B = 0.1 * total_reads_A / total_reads_A  # Here total_reads_A is used to match the 10% of the Sample A size

# Process Sample B: equivalent to 10% of Sample A reads
process_sample("SampleB", fastq1_sampleB_read1, fastq2_sampleB_read2, output_subsetB_read1, output_subsetB_read2, proportion_B)
