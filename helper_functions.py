"""
A file containing some helper functions.
"""

from Bio import SeqIO


def is_fasta(file):
    """
    Check whether a given file is a valid FASTA file.
    https://stackoverflow.com/questions/44293407/how-can-i-check-whether-a-given-file-is-fasta
    """
    with open(file) as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)
