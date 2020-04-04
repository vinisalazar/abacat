#!/usr/bin/env python
"""
Create ANI dendrograms using FastANI and SciPy.

Input is 2 or more gene or genome files.

Example usage as script:

python abacat/dendrogram.py -i my_genomes/

"""

import argparse
import subprocess
import pandas as pd
from os import path, listdir, mkdir
from matplotlib import pyplot as plt
from abacat import prodigal, CONFIG, timer_wrapper
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram, linkage


class ANIDendrogram:
    def __init__(
        self,
        fastani_input="fastani_input.txt",
        threads=CONFIG["threads"],
        cmd=None,
        fraglen=300,
        fastani_output=None,
        output_dir="ani_output/",
        ani_table=None,
        fig_output=None,
        key_df=None,
    ):
        """
        A class to make ANI dendrograms using FastANI and SciPy.

        :param fastani_input: File with path of gene files, one per line, to give to fastANI
        :param threads: Number of threads to use with FastANI
        :param cmd: Command string
        :param fraglen: Length of fragments for FastANI. Default 300
        :param fastani_output: FastANI file name
        :param output_dir: Output directory name
        :param ani_table: Table with FastANI results
        :param fig_output: Dendrogram figure file name
        :param color_threshold: Parameter for color threshold in Dendrogram. Default is 5 for species level.
        :param key_df: Dataframe with two columns: old name and new name, to filter or rename samples in dendrogram.
        """
        super(ANIDendrogram, self).__init__()
        self.fastani_input = fastani_input
        self.threads = threads
        self.cmd = cmd
        self.fraglen = fraglen
        self.fastani_output = fastani_output
        self.output_dir = output_dir
        self.ani_table = ani_table
        self.fig_output = fig_output
        self.fastani_bin = CONFIG["third_party"]["fastANI"]
        self.df = None  # Pandas dataframe which will be saved to ani_table.
        self.key_df = key_df

    def import_key_file(self, key_file):
        self.key_df = import_and_validate_key_file(key_file)

    def make_fastani_input(self, list_of_genomes=None, kind="files"):
        """
        :param list_of_genomes: List of genomes to build file.
        :param kind: 'instances' if instances of Genome class. 'files' if list of files.
        :param run: runs Prodigal for genome if genes file present.
        :return: FastANI input with gene files, one per line.
        """

        if not path.isdir(self.output_dir):
            mkdir(self.output_dir)

        self.fastani_input = path.basename(self.fastani_input)  # Prevents from doubling the directory if running again.
        self.fastani_input = path.join(self.output_dir, self.fastani_input)

        if kind == "instances":
            gene_files = []
            for genome in list_of_genomes:
                if str(type(genome)) != "<class 'abacat.genome.Genome'>":
                    print(f"{genome} is not a valid Genome instance. Ignoring it.")
                    # raise Exception("Please pass a valid Genome instance")
                    pass
                elif not genome.files["prodigal"]["genes"]:
                    raise Exception(f"{genome} does not have a genes file. Ignoring it.")
                else:
                    gene_file = genome.files["prodigal"]["genes"]
                    gene_files.append(gene_file)
            with open(self.fastani_input, "w") as f:
                f.write("\n".join(gene_files))
        elif kind == "files":
            with open(self.fastani_input, "w") as f:
                for genome in list_of_genomes:
                    f"Writing genome to {self.fastani_input}."
                    f.write(genome + "\n")

        if path.isfile(self.fastani_input):
            with open(self.fastani_input) as f:
                for i, l in enumerate(f):
                    pass

            return f"Wrote {i + 1}/{len(list_of_genomes)} gene files to {self.fastani_input}."
        else:
            return f"Could write to {self.fastani_input}."

    def run(self):
        if not self.fastani_output:
            self.fastani_output = path.join(
                self.output_dir, f"fastani_out_{self.fraglen}"
            )

        if not self.cmd:
            self.cmd = f"{self.fastani_bin} --ql {self.fastani_input} --rl {self.fastani_input} -t {self.threads} -o {self.fastani_output} --minFraction 0"

        if self.fraglen:
            self.cmd += f" --fragLen {self.fraglen}"

        @timer_wrapper
        def run_():
            print("Running FastANI.")
            subprocess.call(self.cmd, shell=True)

        run_()
        if path.isfile(self.fastani_output):
            print(f"FastANI ran successfully!")

    def make_ani_table(self):
        """
        :return: ANI distance table from FastANI output.
        """

        # Setting the columns and reading the dataframe
        columns = [
            "Genome_A",
            "Genome_B",
            "ANI",
            "orthologous_fraction",
            "total_fragments",
        ]
        df = pd.read_csv(self.fastani_output, sep="\t", header=None)
        df.columns = columns

        # Creating the concatenated data frame
        df_ = df.copy()
        df_.rename(columns={"Genome_A": "Genome_B", "Genome_B": "Genome_A"}, inplace=True)
        df_t = pd.concat([df, df_], sort=True)
        df_t.ANI = round(df_t.ANI, 2)

        # Switching from 'long' form to a distance matrix
        df = pd.pivot_table(
            df_t, values="ANI", index=("Genome_B"), columns=("Genome_A")
        )
        del df_
        del df_t  # Don't need them anymore

        # Make sure the horizontal diagonal is symmetrical
        for i, j in zip(range(len(df)), range(len(df))):
            df.iloc[i, j] = round(99.99, 2)

        df.dropna(inplace=True)  # Drop na columns to prevent

        """
        The df is using the full path and genes file.
        This next bit replaces it with a shorter, prettier name.
        And allows to filter and rename with the key file.
        """

        df.index = [path.basename(i).replace("_prodigal_genes", "") for i in df.index]
        df.columns = [
            path.basename(i).replace("_prodigal_genes", "") for i in df.columns
        ]

        if not self.ani_table:
            self.ani_table = path.join(self.output_dir, "ani_table.txt")

        self.df = df
        print(
            f"Writing ANI distance matrix to {path.abspath(self.ani_table)}. It is {df.shape[0]} by {df.shape[1]}."
        )
        df.to_csv(self.ani_table, header=False, index=False)

    def make_dendrogram(self, color_threshold=5, filter_rename=False):
        """
        :param filter_rename: Will filter and rename using self.key_file. If there is not a key_file, will attempt to import from passed key_file.
        :param color_threshold: Color threshold to paint dendrogram branches.
        :return: Build dendrogram.
        """
        print(f"Plotting dendrogram with color threshold of {color_threshold}.")
        if not any(self.df):
            raise Exception(
                "You don't have an ANI table to make a dendrogram. Run make_ani_table() first."
            )

        table = self.df
        if filter_rename:
            if not path.isfile(filter_rename):
                self.import_key_file(filter_rename)

            table.rename(
                columns=dict(zip(self.key_df["old_name"], self.key_df["new_name"])),
                index=dict(zip(self.key_df["old_name"], self.key_df["new_name"])),
                inplace=True,
            )
            table = table.loc[
                list(self.key_df["new_name"]), list(self.key_df["new_name"])
            ]
        table = abs(table - 99.99)
        X = squareform(table)
        Z = linkage(X, method="complete", metric="cityblock", optimal_ordering=True)

        fig, ax = plt.subplots()

        augmented_dendrogram(
            Z,
            labels=table.columns,
            leaf_rotation=-90,
            color_threshold=color_threshold,
            leaf_font_size=12,
            ax=ax,
        )

        if not self.fig_output:
            self.fig_output = path.join(
                self.output_dir, f"ANI_dendrogram_ct_{color_threshold}"
            )
        plt.savefig(self.fig_output, bbox_inches="tight")
        if path.isfile(self.fig_output + ".png"):
            print(f"Generated dendrogram at {path.abspath(self.fig_output)}.")


def augmented_dendrogram(*args, **kwargs):

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get("no_plot", False):
        for i, d in zip(ddata["icoord"], ddata["dcoord"]):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > 1.5:
                plt.plot(x, y, "ro")
                plt.annotate(
                    "%.3g" % y,
                    (x, y),
                    xytext=(0, -8),
                    textcoords="offset points",
                    va="top",
                    ha="center",
                )

    return ddata


def import_and_validate_key_file(key_file):
    """
    :param key_file: Validates key file to make sure it will work.
    :return: Pandas dataframe with key file
    """
    with open(key_file) as f_:
        r = [i.split() for i in f_]

    for i in r:
        if len(i) != 2:
            raise Exception(
                "Invalid key file. Must be a tab delimited file with two columns, old name and new name."
            )

    r = pd.DataFrame(r)
    r.columns = ("old_name", "new_name")

    return r


if __name__ == "__main__":

    @timer_wrapper
    def main():
        parser = argparse.ArgumentParser(
            description="Plot a dendrogram based on ANI distances."
        )
        parser.add_argument(
            "-i",
            "--input",
            help="Input genomes. "
            "May be a file containing the path to each genome, one per line, or a folder containing genomes. "
            "Make sure to flag genes=True if using file containing gene sequences "
            "(otherwise will attempt to run Prodigal before).",
        )
        parser.add_argument(
            "-o",
            "--output",
            help="Path to the output directory. Will create 4 files:\n"
            "1. FastANI input file (skipped if input is file containing one genome per line).\n"
            "2. FastANI output file (tab-delimited).\n"
            "3. ANI distance matrix.\n"
            "4. Dendrogram plot.\n",
            default="./fastani_output",
        )
        parser.add_argument(
            "-g",
            "--genes",
            help="Flag to assign whether we are providing whole genome (False) or only genes (True). "
            "If False, will run Prodigal to predict genes before running FastANI (default: False).",
            type=bool,
            default=False,
        )
        parser.add_argument(
            "-k",
            "--keys",
            help="Use this parameter to rename or filter your samples. "
            "Must be a space or tab delimited file of <OLD-NAME> <NEW NAME>."
            "Samples not in this file will be filtered out.",
            default=None,
        )
        parser.add_argument(
            "-ct",
            "--color_threshold",
            help="Color threshold to paint branches of dendrogram. Default is 5 (species level).",
            type=int,
            default=5,
        )
        parser.add_argument(
            "-fl",
            "--fragLen",
            help="Fragment length to give to FastANI. Default is 300.",
            type=int,
            default=300,
        )
        parser.add_argument(
            "-t",
            "--threads",
            help="Number of threads to use. Default is set in Abacat config file.",
            type=int,
            default=CONFIG["threads"],
        )
        parser.add_argument(
            "-s",
            "--skip",
            help="Skip FastANI run. Only plot figure (must provide FastANI output).",
            type=bool,
            default=False,
        )

        args = parser.parse_args()

        # Parsing args.input
        try:
            if path.isdir(args.input):
                input_files = [
                    path.abspath(path.join(args.input, i)) for i in listdir(args.input)
                ]
                input_files = [
                    i for i in input_files if i.endswith((".fna", ".fasta", ".fa"))
                ]  # Only catch some extensions.
                if args.genes:
                    input_files = [i for i in input_files if "prodigal_genes" in i]
                else:
                    input_files = [
                        i for i in input_files if "prodigal" not in i
                    ]  # Skip previous Prodigal files
            elif path.isfile(args.input):
                if not args.genes:
                    with open(args.input) as f:
                        input_files = [path.abspath(i.strip()) for i in f.readlines()]
                else:
                    input_files = (
                        args.input
                    )  # If it is already a genes file, it is proper for FastANI. No need to write another.
        except FileNotFoundError as error:
            raise (error)
            input_files = None  # This line is only to suppress warnings.

        # Parsing args.ouput
        print(f"Set {path.abspath(args.output)} as output directory.")
        if not path.isdir(args.output):
            mkdir(args.output)
        fastani_input = path.join(args.output, "fastani_input.txt")

        # Parsing args.genes -- Prodigal runs
        if not args.genes:
            print(f"You have {len(input_files)} genomes to be processed:\n")
            print("\n".join((path.basename(i) for i in input_files)), "\n")
            print(f"Starting gene prediction.\n")
            ix, gene_files = 1, []
            for file in input_files:
                print(f"Processing file {ix}/{len(input_files)}.")
                genes = prodigal.run(file, output=path.dirname(file), quiet=True, print_files=False)
                gene_files.append(genes["genes"])
                ix += 1
            with open(fastani_input, "w") as f:
                f.write("\n".join(gene_files))
        else:
            with open(fastani_input, "w") as f:
                f.write("\n".join(input_files))

        ani = ANIDendrogram(
            output_dir=args.output,
            fastani_input=fastani_input,
            fastani_output=path.join(
                args.output, f"fastani_out_{args.fragLen}.tsv"
            ),
            threads=args.threads,
        )
        if not args.skip:
            ani.run()
        ani.make_ani_table()
        if args.keys:
            print(f"Key file set as {path.abspath(args.keys)}.")
            ani.make_dendrogram(
                color_threshold=args.color_threshold, filter_rename=args.keys
            )
        else:
            ani.make_dendrogram(color_threshold=args.color_threshold)

    # https://www.youtube.com/watch?v=odeHP8N4LKc
    main()
