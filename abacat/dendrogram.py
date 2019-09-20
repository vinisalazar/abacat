"""
Create ANI dendrogram using FastANI and Scipy.
"""

import argparse
import subprocess
import pandas as pd
from os import path, listdir
from matplotlib import pyplot as plt
from abacat import CONFIG, timer_wrapper
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram, linkage


class ANI_dendrogram:
    def __init__(
        self,
        fastani_input="fastani_input.txt",
        threads=CONFIG["threads"],
        cmd=None,
        fraglen=3000,
        minfrag=50,
        fastani_output=None,
        fastani_dir=CONFIG["data_dir"],
        ani_table=None,
        fig_output=None,
        color_threshold=30
    ):
        """
        :param fastani_input: File with path of gene files, one per line, to give to fastANI
        :param threads: Number of threads to use with FastANI
        :param cmd: Command string
        :param fraglen: Length of fragments for FastANI. Default 3000
        :param minfrag: Minimum number of fragments to trust ANI. Default 50
        :param fastani_output: Output file name
        :param ani_table: Table with FastANI results
        :param fig_output: Dendrogram figure file name
        :param color_threshold: Parameter for color threshold in Dendrogram. Default is 30 for genus level.
        """
        super(ANI_dendrogram, self).__init__()
        self.fastani_input = fastani_input
        self.threads = threads
        self.cmd = cmd
        self.fraglen = fraglen
        self.minfrag = minfrag
        self.fastani_output = fastani_output
        self.fastani_dir = fastani_dir
        self.ani_table = ani_table
        self.fig_output = fig_output
        self.color_threshold = color_threshold
        self.df = None  # Pandas dataframe which will be saved to ani_table.


    def make_fastani_input(self, list_of_genomes=None, kind="instances"):
        """
        :param output_dir: Directory to create list of files.
        :param list_of_genomes: List of genomes to build file.
        :kind: 'instances' if instances of Genome class. 'files' if list of files.
        :return: FastANI input with gene files, one per line.
        """

        filename = path.join(self.fastani_dir, self.fastani_input)

        if kind == "instances":
            for genome in list_of_genomes:
                if str(type(genome)) != "<class 'abacat.genome.Genome'>":
                    print(f"{genome} is not a valid Genome instance. Ignoring it.")
                    #raise Exception("Please pass a valid Genome instance")
                    pass
                elif not genome.files["prodigal"]["genes"]:
                    print(f"{genome} does not have a genes file. Ignoring it.")
                    pass
                else:
                    gene_file = genome["prodigal"]["genes"]
                    with open(filename, "w") as f:
                        f.write(gene_file + "\n")
        elif kind == "files":
            for genome in list_of_genomes:
                with open(filename) as f:
                    f.write(genome + "\n")

        if path.isfile(filename):
            with open(filename) as f:
                for i, l in enumerate(f):
                    pass

                return f"Wrote {i}/{len(list_of_genomes)} gene files to {filename}."

    def run(self):
        if not self.fastani_output:
            self.fastani_output = path.join(self.fastani_dir, f"fastani_out_{self.fraglen}_{self.fraglen}")

        if not self.cmd:
            self.cmd = f"FastANI --ql {self.fastani_input} --rl {self.fastani_input} -t {self.threads} -o {self.fastani_output}"

        if self.fraglen:
            self.cmd += f" --fraglen {self.fraglen}"
        if self.minfrag:
            self.cmd += f" --minfrag {self.minfrag}"

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
        columns = ["Genome_A", "Genome_B", "ANI", "orthologous_fraction", "total_fragments"]
        df = pd.read_csv(self.fastani_output, columns=columns)
        df_ = df
        df_.rename(columns={"Genome_A": "Genome_B", "Genome_B": "Genome_A"}, inplace=True)
        df = pd.concat([df, df_], sort=True)
        del df_
        df.ANI = round(df.ANI, 2)
        df = pd.pivot_table(df, values="ANI", index=("Genome_A"), columns=("Genome_B"))

        for i, j in zip(range(len(df)), range(len(df))):
            df.iloc[i, j] = round(99.99, 2)

        if not self.ani_table:
            self.ani_table = "ani_table.txt"

        self.ani_table = path.join(self.fastani_dir, self.ani_table)
        self.df = df
        print(f"Writing ANI table to {self.ani_table}. It is {df.shape[0]} by {df.shape[1]}.")
        df.to_csv(self.ani_table, header=False, index=False)

    def make_dendrogram(self):
        """
        :return: Build dendrogram.
        """
        X = squareform(abs(self.df - 99.99))
        Z = linkage(X, method="complete", metric="cityblock", optimal_ordering=True)

        fig, ax = plt.subplots()

        augmented_dendrogram(Z, labels = self.df.columns, leaf_rotation=-90, color_threshold = self.color_threshold, leaf_font_size=12, ax=ax)

        if not self.fig_output:
            self.fig_output(f"ANI_dendrogram_ct_{self.color_threshold}")

        self.fig_output = path.join(self.fastani_dir, self.fig_output)
        plt.savefig(self.fig_output, bbox_inches="tight")


def augmented_dendrogram(*args, **kwargs):

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        for i, d in zip(ddata['icoord'], ddata['dcoord']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > 1.5:
                plt.plot(x, y, 'ro')
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -8),
                         textcoords='offset points',
                         va='top', ha='center')

    return ddata


if __name__ == "__main__":
    run_dir = "/Users/viniWS/Bio/abacat_data/synecho/fastani/"
    list_of_genomes = [path.join(run_dir, i) for i in run_dir]
    ani = ANI_dendrogram()
    ani.make_fastani_input(kind="files", list_of_genomes=list_of_genomes)
    ani.run()
    ani.make_dendrogram()