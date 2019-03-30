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
import argparse
import subprocess


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


def prodigal(file, type="np", output=None, directory=None):
    """
    Calls Prodigal on an input file.

    Input:
    .fasta or .fna file.
    Output:
    Genes (.fna) and proteins (.faa) files.
    """

    if not output:
        output = os.path.basename(file)
        output = os.path.splitext(output)[0]

    prodigal = subprocess.Popen(
        f"prodigal -i {file} -a {file_name + '_proteins.faa'} \
        -d {file_name + '_genes.fna'} -o {file_name + '_out.txt'} \
        -s {file_name + '_scores.txt'}", shell=True, stdout=subprocess.PIPE
    )

    # TRY TO FIX THIS LATER
    # redirection tip came from here
    # https://eli.thegreenplace.net/2015/redirecting-all-kinds-of-stdout-in-python/
    # prodigal_output = prodigal.communicate()[0]
    #
    # with open(f"{file_name}_prodigal.log", "w") as f:
    #     f.write(str(prodigal_output))
