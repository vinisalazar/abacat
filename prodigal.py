"""
A script to call Prodigal to predict both genes and proteins.

Input:
File or folder with contigs in fasta or fna format.

Output:
Separate files or folders with genes and proteins (default) or only one of each.
Includes log and stats file (errors, exceptions, how many genes and proteins)
When output is multiple, includes single log file.

Example usage:

    python prodigal.py -i contigs/
    python prodigal.py -i contigs.fasta -o proteins.faa -t n -l log_stats.txt
"""

import os
import sys
import argparse
import subprocess
from Bio import SeqIO


def get_input(input_path):
    """
    Checks if the input is a file or directory.
    Returns the file path and type.
    """
    if os.path.isfile(input_path):
        return os.path.abspath(input_path), "file"

    elif os.path.isdir(input_path):
        return os.path.abspath(input_path), "dir"

    else:
        # Change this to return None?
        raise FileNotFoundError


def prodigal(file, output=None):
    """
    Calls Prodigal on an input file.

    Input
    .fasta or .fna file.
    Output
    Genes (.fna), proteins (.faa), gene scores (.txt), gbk file (.gbk).
    """

    # Check if it is a valid fasta file.
    # https://stackoverflow.com/questions/44293407/how-can-i-check-whether-a-given-file-is-fasta
    def is_fasta(file):
        with open(file) as handle:
            fasta = SeqIO.parse(handle, "fasta")
            return any(fasta)

    if not is_fasta(file):
        raise Exception(
            "Your file is not valid. Please check if it is a valid FASTA file."
        )

    # Create default output format
    if not output:
        output = os.path.splitext(file)[0]
    else:
        output = output

    if not os.path.isdir(output):
        os.mkdir(output)

    output = os.path.join(output, output.split('/')[-1])

    prodigal = subprocess.call(
        f"prodigal -i {file} -a {output + '_proteins.faa'} \
        -d {output + '_genes.fna'} -o {output + '_cds.gbk'} \
        -s {output + '_scores.txt'}", shell=True
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""
    A script to call Prodigal to predict both genes and proteins.\n

    Input:\n
    File or folder with contigs in fasta or fna format.\n

    Output:\n
    Separate files or folders with genes and proteins (default) \n
    Includes log and stats file (errors, exceptions, how many genes and proteins)\n
    When output is multiple, includes single log file.\n

    Example usage:\n

        python prodigal.py -i contigs/\n
        python prodigal.py -i contigs.fasta -o proteins.faa\n
    """
    )

    parser.add_argument("-i", "--input", help="Input file. Must be a valid\
                        FASTA file.")
    parser.add_argument("-o", "--output", help="Path to output folder",
                        default="")

    args = parser.parse_args()

    if not args.input:
        parser.print_help()
        sys.exit(0)

    prodigal(args.input, args.output)
