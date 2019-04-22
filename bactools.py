"""
A file containing our main classes and functions.
"""

import datetime
import time
from Bio import SeqIO


class Assembly(object):
    """
    Assembly, a class containing bacterial assembly data and methods.
    """

    def __init__(self, arg):
        super(Assembly, self).__init__()
        self.arg = arg


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
