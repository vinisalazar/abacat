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


def is_fasta_wrapper(func):
    """
    A decorator wrap the functions which need checking for a valid FASTA.
    """
    def wrapper(*args, **kwargs):
        if not is_fasta(args[0]):
            raise Exception(
                "Your file is not valid. Please check if it is a valid FASTA file."
            )
        func(args(*args, **kwargs))
    return wrapper
