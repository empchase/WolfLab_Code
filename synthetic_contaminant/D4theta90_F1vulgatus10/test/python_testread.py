import gzip

def read_gzip_file(file_path):
    with gzip.open(file_path, 'rt') as f:
        for _ in range(5):
            print(f.readline())

read_gzip_file('/global/scratch/projects/fc_wolflab/bbioinfo/empchase/synthetic_contaminant/subsetD4_read1.fastq.gz')
read_gzip_file('/global/scratch/projects/fc_wolflab/bbioinfo/empchase/synthetic_contaminant/subsetD4_read2.fastq.gz')
