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
    python prodigal.py -i contigs.fasta
"""

import os
import sys
import time
import argparse
import datetime
import subprocess
from bactools_helper import is_fasta_wrapper, timer_wrapper


@is_fasta_wrapper
def prodigal(file, output=None, quiet=False):
    """
    Calls Prodigal on an input file.

    Input:
    Valid FASTA file.
    Output:
    Genes (.fna), proteins (.faa), gene scores (.txt), gbk file (.gbk).

    Example:

        from prodigal import prodigal
        prodigal("contigs.fasta", "output_folder/")

    """

    if not output:
        output = os.path.splitext(file)[0] + "_prodigal"
    else:
        output = os.path.join(
            output, os.path.basename(os.path.splitext(file)[0]) + "_prodigal"
        )

    if not os.path.isdir(output):
        os.mkdir(output)

    output = os.path.join(os.path.abspath(output), output.split("/")[-1])
    output_files = {
        "genes": output + "_genes.fna",
        "proteins": output + "_proteins.faa",
        "cds": output + "_cds.gbk",
        "scores": output + "_scores.txt",
    }

    cmd = f"prodigal -i {file} -a {output_files['proteins']} \
            -d {output_files['genes']} -o {output_files['cds']} \
            -s {output_files['scores']}"

    # This suppresses console output from Prodigal
    if quiet:
        cmd += " -q"

    prodigal = subprocess.call(cmd, shell=True)
    print(f"Created files at {output}:")
    for k, v in output_files.items():
        print("\t", v)

    return output_files


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="""
    A script to call Prodigal to predict both genes and proteins.
    """
    )

    parser.add_argument(
        "-i", "--input", help="Input FASTA file or dir containing fasta files"
    )
    parser.add_argument("-o", "--output", help="Path to output folder", default=".")

    args = parser.parse_args()

    if not args.input:
        parser.print_help()
        sys.exit(0)

    input = args.input

    @timer_wrapper
    def main():
        if os.path.isfile(input):
            print(f"Starting script. Your input file is {input}.")
            prodigal(input, args.output)

        elif os.path.isdir(input):

            files = os.listdir(input)
            files = [os.path.join(input, i) for i in files]
            files = [i for i in files if os.path.isfile(i)]

            print("\n")
            print(
                f"Starting Prodigal. Quiet setting is {quiet}. You have {len(files)} files to be processed in {input}:\n"
            )
            print("\n".join(files), "\n")

            success = 0
            failure = 0

            for i in files:
                try:
                    print(f"Running Prodigal for {i}.")
                    prodigal(i, args.output, quiet=True)
                    if os.path.isdir(os.path.splitext(i)[0]):
                        success += 1
                except Exception as err:
                    print(f"Error for {i}. Please check if it is a valid FASTA file.")
                    failure += 1
                    pass

            print("\n")
            print(f"Done. {success} assemblies processed. {failure} errors.")

        else:
            raise FileNotFoundError

    main()
