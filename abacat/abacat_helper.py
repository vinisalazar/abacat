from Bio import SeqIO
import datetime
import time
import logging


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
        try:
            if not is_fasta(args[0]):
                raise Exception(
                    "Your file is not valid. Please check if it is a valid FASTA file."
                )
            return func(*args, **kwargs)
        except IndexError:
            logging.error("Please specify an input file.", exc_info=True)

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
        logging.info(f"Took {delta}")

    return wrapper


@is_fasta_wrapper
def get_records(fasta_file, kind="gen"):
    with open(fasta_file) as f:
        if kind == "gen":
            records = SeqIO.parse(fasta_file, format="fasta")
        elif kind == "list":
            records = list(SeqIO.parse(fasta_file, format="fasta"))
        elif kind == "dict":
            records = SeqIO.to_dict(SeqIO.parse(fasta_file, format="fasta"))
        else:
            logging.info(f"Specified {kind} kind. Please specify a valid kind.")

    return records
