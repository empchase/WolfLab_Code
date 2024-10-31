# test_biopython.py
try:
    from Bio import SeqIO
    print("Biopython is installed and working.")
except ImportError:
    print("Biopython is not installed.")
