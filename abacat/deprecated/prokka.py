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
from abacat.abacat_helper import is_fasta_wrapper, timer_wrapper


@is_fasta_wrapper
def prokka(fasta_file, output="", cpus=None, params=""):
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

    if not fasta_file:
        return "Please specify an input file."

    if output:
        output = os.path.join(os.path.abspath(output), os.path.splitext(fasta_file)[0])
    else:
        output = os.path.basename(os.path.splitext(fasta_file)[0])
    prefix = f"--prefix {output}"
    output = f"{output}_prokka"
    output_ = f"--outdir {output}"

    if not cpus:
        cpus = int(os.cpu_count() / 2)
    cpus = f"--cpus {cpus}"

    cmd = f"prokka {output_} {prefix} {fasta_file} {cpus} {params}"

    prokka = subprocess.call(cmd, shell=True)

    output_files = dict()

    for f in os.listdir(output):
        f = os.path.join(os.path.abspath(output), f)
        ext = os.path.splitext(f)[1]
        # How to change this to a more elegant dictionary?
        if ext == ".fna":
            output_files["input_contigs"] = f
        elif ext == ".faa":
            output_files["proteins"] = f
        elif ext == ".ffn":
            output_files["genes"] = f
        elif ext == ".fsa":
            output_files["submit_contigs"] = f
        elif ext == ".tbl":
            output_files["featuretable"] = f
        elif ext == ".sqn":
            output_files["sequin"] = f
        elif ext == ".gbk":
            output_files["genbank"] = f
        elif ext == ".gff":
            output_files["gff"] = f
        elif ext == ".log":
            output_files["log"] = f
        elif ext == ".txt":
            output_files["stats"] = f

    return output_files


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="""
    A script to call Prokka for genome annotation.
    """
    )

    parser.add_argument(
        "-i", "--input", help="Input FASTA file or dir containing fasta files."
    )
    parser.add_argument("-o", "--output", help="Path to output folder.", default="")
    parser.add_argument(
        "-p",
        "--params",
        help="Additional parameters to pass to Prokka as a string.",
        default="",
    )
    parser.add_argument(
        "-c",
        "--cpus",
        help="Number of CPUs to use. Default is half of the system's CPUs.",
        default=None,
    )

    args = parser.parse_args()

    if not args.input:
        parser.print_help()
        sys.exit(0)

    input = args.input

    @timer_wrapper
    def main():
        if os.path.isfile(input):
            print(f"Starting script. Your input file is {input}.")
            prokka(input, args.output, args.cpus, args.params)

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
                    prokka(i, args.output, args.cpus, args.params)
                    if os.path.isdir(os.path.splitext(i)[0]):
                        success += 1
                except Exception as err:
                    print(f"Error for {i}. Please check if it is a valid FASTA file.")
                    failure += 1
                    pass

        else:
            raise FileNotFoundError

    main()
