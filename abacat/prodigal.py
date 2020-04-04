!#/usr/bin/env python
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
import argparse
import subprocess
from abacat.abacat_helper import is_fasta, is_fasta_wrapper, timer_wrapper


class Prodigal:
    """
    This class will hold Prodigal run data.

    Calls Prodigal on an input file.

    Input:
    Valid FASTA file.
    Output:
    Genes (.fna), proteins (.faa), gene scores (.txt), gbk file (.gbk).

    """

    def __init__(self, contigs, output=None, quiet=False, scores=False):
        super(Prodigal, self).__init__()
        self.name = None
        is_fasta(contigs)
        self.contigs = contigs  # an assembled genome contigs file.
        self.quiet = quiet
        self.finished = None
        self.scores = scores

        if not output:
            output = os.path.join(
                os.getcwd(),
                os.path.basename(os.path.splitext(self.contigs)[0]) + "_prodigal",
            )
        else:
            output = os.path.join(
                output,
                os.path.basename(os.path.splitext(self.contigs)[0]) + "_prodigal",
            )
        self.output = os.path.join(os.path.abspath(output), output.split("/")[-1])
        self.output_files = {
            "genes": output + "_genes.fna",
            "proteins": output + "_proteins.faa",
            "cds": output + "_cds.gbk",
        }

        self.cmd = f"prodigal -i {self.contigs} -a {self.output_files['proteins']} \
                    -d {self.output_files['genes']} -o {self.output_files['cds']}"

        if self.scores:
            self.output_files["scores"] = output + "_scores.txt"
            self.cmd += f" -s {self.output_files['scores']}"
        if self.quiet:
            self.cmd += " -q"

    def run(self, print_files=False):
        is_fasta_wrapper(self.contigs)
        subprocess.call(self.cmd, shell=True)

        if all(os.path.isfile(value) for _, value in self.output_files.items()):
            self.finished = True
            if print_files:
                print(f"Created files at {self.output}:")
                for _, v in self.output_files.items():
                    print("\t", v)

        else:
            self.finished = False

        return self.output_files


def run(contig_file, output=None, quiet=False, print_files=False):
    """
    Run outside of class scope.
    """
    p = Prodigal(contig_file, output=output, quiet=quiet)
    p.run(print_files=print_files)

    return p.output_files


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

    input_ = args.input

    @timer_wrapper
    def main():
        if os.path.isfile(input_):
            print(f"Starting script. Your input file is {input}.")
            p = Prodigal(input_, output=args.output)
            p.run()

        elif os.path.isdir(input_):

            files = os.listdir(input_)
            files = [os.path.join(input_, i) for i in files]
            files = [i for i in files if os.path.isfile(i)]

            print("\n")
            print(
                f"Starting Prodigal. You have {len(files)} files to be processed in {input_}:\n"
            )
            print("\n".join(files), "\n")

            success = 0
            failure = 0

            for contig_file in files:
                try:
                    print(f"Running Prodigal for {contig_file}.")
                    run(contig_file)
                    if os.path.isdir(os.path.splitext(contig_file)[0]):
                        success += 1
                except ValueError:
                    print(
                        f"Error for {contig_file}. Please check if it is a valid FASTA contigs file."
                    )
                    failure += 1
                    pass

            print("\n")
            print(f"Done. {success} assemblies processed. {failure} errors.")

        else:
            raise FileNotFoundError

    main()
