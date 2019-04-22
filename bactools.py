"""
A file containing our main classes and functions.
"""

import datetime
import os
import time
from Bio import SeqIO
from prodigal import prodigal


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
    A decorator wrapper for functions which need checking for a valid FASTA.
    """

    def wrapper(*args, **kwargs):
        if not is_fasta(args[0]):
            raise Exception(
                "Your file is not valid. Please check if it is a valid FASTA file."
            )
        return func(*args, **kwargs)

    return wrapper


def timer_wrapper(func):
    """
    A wrapper to time the execution time of our functions.

    Example usage:

    if __name__ == "__main__":
        *define argparse here*

        @timer_wrapper
        def main():
            '''
            This is your main function
            '''

        main()
    """

    def wrapper(*args, **kwargs):
        start = time.time()

        func(*args, **kwargs)

        end = time.time()
        delta = str(datetime.timedelta(seconds=end - start))
        print(f"Took {delta}")

    return wrapper


class Assembly:
    """
    Assembly, a class containing bacterial assembly data and methods.
    """

    def __init__(self):
        super(Assembly, self).__init__()
        self.inputs = dict()
        self.outputs = dict()
        self.metadata = None
        self.geneset = None
        self.protset = None

    @timer_wrapper
    def run_prodigal(self, quiet=True):
        if not self.inputs["contigs"]:
            raise Exception("Must specify input contigs file.")

        print(f"Starting script. Your input file is {input}.")

        prodigal_out = prodigal(self.inputs["contigs"], quiet=quiet)

        self.outputs["prodigal"] = prodigal_out

    def load_from_fasta(self, input_fasta, metadata=None):
        self.inputs["contigs"] = os.path.abspath(input_fasta)
