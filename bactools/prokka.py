"""
A script to call Prokka to do it's things. Similar to prodigal.py.

Input:
file or folder with contigs in fasta or fna format.

Output:
Prokka folder to selected output (default current directory).

Example usage:

    python prokka.py -i contigs/
    python prokka.py -i contigs.fasta

"""

import os
import sys
import time
import argparse
import datetime
import subprocess
from .bactools_helper import is_fasta_wrapper, timer_wrapper


@is_fasta_wrapper
def prokka(file, output="", cpus=None):
    """
    Calls Prokka on an input file.

    Input:
    Valid FASTA file.

    Output:
    Prokka output to selected folder.

    Example:

        from prokka import prokka
        prokka("contigs.fasta", "output_folder/")
    """

    if output:
        output = os.path.join(os.path.abspath(output), os.path.splitext(file)[0])
    else:
        output = os.path.basename(os.path.splitext(file)[0])
    output = f"--outdir {output}_prokka"

    if not cpus:
        cpus = int(os.cpu_count() / 2)
    cpus = f"--cpus {cpus}"

    cmd = f"prokka {output} {file} {cpus}"

    prokka = subprocess.call(cmd, shell=True)

    # Build dict for return value.
    return 0


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="""
    A script to call Prokka for genome annotation.
    """
    )

    parser.add_argument(
        "-i", "--input", help="Input FASTA file or dir containing fasta files"
    )
    parser.add_argument("-o", "--output", help="Path to output folder", default="")

    args = parser.parse_args()

    if not args.input:
        parser.print_help()
        sys.exit(0)

    input = args.input

    @timer_wrapper
    def main():
        if os.path.isfile(input):
            print(f"Starting script. Your input file is {input}.")
            prokka(input, args.output)

        elif os.path.isdir(input):

            files = os.listdir(input)
            files = [os.path.join(input, i) for i in files]
            files = [i for i in files if os.path.isfile(i)]

            print("\n")
            print(
                f"Starting script. You have {len(files)} files to be processed in {input}:\n"
            )
            print("\n".join(files), "\n")

            success = 0
            failure = 0

            for i in files:
                try:
                    print(f"Running Prokka for {i}.")
                    prokka(i, args.output)
                    if os.path.isdir(os.path.splitext(i)[0]):
                        success += 1
                except Exception as err:
                    print(f"Error for {i}. Please check if it is a valid FASTA file.")
                    failure += 1
                    pass

        else:
            raise FileNotFoundError

    main()
