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
from bactools.bactools_helper import is_fasta_wrapper, timer_wrapper


class Prodigal:
    """
    This class will hold Prodigal run data.

    Calls Prodigal on an input file.

    Input:
    Valid FASTA file.
    Output:
    Genes (.fna), proteins (.faa), gene scores (.txt), gbk file (.gbk).

    Example:

    from prodigal import prodigal
    prodigal("contigs.fasta", "output_folder/")

    """

    def __init__(self, contigs, output=None, quiet=False):
        super(Prodigal, self).__init__()
        self.name = None
        self.contigs = contigs  # an assembled genome contigs file.
        self.quiet = quiet
        self.finished = None

        if not output:
            output = os.path.join(
                os.getcwd(), os.path.basename(os.path.splitext(self.contigs)[0]) + "_prodigal"
            )
        else:
            output = os.path.join(
                output, os.path.basename(os.path.splitext(self.contigs)[0]) + "_prodigal"
            )
        self.output = os.path.join(os.path.abspath(output), output.split("/")[-1])
        self.output_files = {
            "genes": output + "_genes.fna",
            "proteins": output + "_proteins.faa",
            "cds": output + "_cds.gbk",
            "scores": output + "_scores.txt",
        }

        self.cmd = f"prodigal -i {self.contigs} -a {self.output_files['proteins']} \
                    -d {self.output_files['genes']} -o {self.output_files['cds']}"

        if self.quiet:
            self.cmd += " -q"

    def run(self):
        is_fasta_wrapper(self.contigs)
        subprocess.call(self.cmd, shell=True)

        if all(os.path.isfile(value) for _, value in self.output_files.items()):
            self.finished = True
            print(f"Created files at {self.output}:")
            for k, v in self.output_files.items():
                print("\t", v)


        else:
            self.finished = False

        return self.output_files



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
            p = Prodigal(input, output = args.output)
            p.run()

        elif os.path.isdir(input):

            files = os.listdir(input)
            files = [os.path.join(input, i) for i in files]
            files = [i for i in files if os.path.isfile(i)]

            print("\n")
            print(
                f"Starting Prodigal. You have {len(files)} files to be processed in {input}:\n"
            )
            print("\n".join(files), "\n")

            success = 0
            failure = 0

            for contig_file in files:
                try:
                    print(f"Running Prodigal for {i}.")
                    @timer_wrapper
                    def run(contig_file):
                        p = Prodigal(contig_file, output=args.output, quiet=True)
                        p.run()

                    run(contig_file)
                    if os.path.isdir(os.path.splitext(contig_file)[0]):
                        success += 1
                except Exception as err:
                    print(f"Error for {contig_file}. Please check if it is a valid FASTA contigs file.")
                    failure += 1
                    pass

            print("\n")
            print(f"Done. {success} assemblies processed. {failure} errors.")

        else:
            raise FileNotFoundError

    main()
