"""
A file containing our main classes and functions.
# TODO: Create class for sets (geneset, protset)
"""

import os
import json
import logging
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import (
    NcbiblastnCommandline,
    NcbiblastpCommandline,
    NcbiblastxCommandline,
)
from abacat.abacat_helper import get_records, is_fasta, is_fasta_wrapper, timer_wrapper
from abacat.prodigal import Prodigal
from abacat.deprecated import prokka
from abacat.config import CONFIG, pathways

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)

logger = logging.getLogger(__name__)


class Genome:
    """
    Genome, a class containing bacterial genome data and methods.

    # TODO: seqstats or Prodigal basic info (length, GC content) for metadata
    """

    def __init__(self, contigs=None, prodigal=False, name=None, directory=None):
        super(Genome, self).__init__()
        self.name = name
        self.directory = directory
        self.seqstats = None
        self.files = dict()
        self.geneset = dict()
        self.protset = dict()
        self.pathways = None

        if contigs:
            self.load_contigs(contigs)

        if prodigal:
            self.load_prodigal()

    """
    Sequence stats methods.
    """

    def print_seqstats(self):
        if self.seqstats:
            print("Sequence stats:\n")
            for key, value in self.seqstats.items():
                print(f"{key}\t\t\t{value}")
        else:
            try:
                self.load_seqstats()
                self.seqstats
            except:
                logger.error(
                    "Tried loading seqstats, but an error occurred.", exc_info=True
                )
                raise

    def fnoseqs(self, kind="prodigal", seqs="genes"):
        """
        Fast check of number of sequences in file.
        """
        with open(self.files[kind][seqs]) as f:
            records = SeqIO.parse(f, "fasta")
            len_ = sum(1 for x in records)

        return len_

    def sseqs(self, kind="prodigal", seqs="genes"):
        """
        Fast check of size of sequences file.
        from https://stackoverflow.com/questions/12523586/python-format-size-application-converting-b-to-kb-mb-gb-tb/37423778
        """
        size = os.stat(self.files[kind][seqs]).st_size

        def format_bytes(size):
            power = 2 ** 10
            n = 0
            power_labels = {0: "", 1: "k", 2: "m", 3: "g", 4: "t"}
            while size > power:
                size /= power
                n += 1

            return f"{round(size, 2)} {power_labels[n]}B"

        return format_bytes(size)

    def valid_contigs(self, quiet=True):
        try:
            any(self.files["contigs"])
        except Exception as e:
            logger.error("Must specify input contigs file.")
        else:
            if not quiet:
                logger.info("Your contig file is a valid FASTA file.")
            return True

    """
    Methods to load data into class instances.
    """

    def load_contigs(self, contigs):
        """
        This loads the contigs input file. I want to make it a constructor method.
        """
        contigs = os.path.abspath(contigs)

        # Must use class decorator, same as in Prodigal class.
        if not is_fasta(os.path.abspath(contigs)):
            raise Exception(
                "Your file is not valid. Please check if it is a valid FASTA file."
            )
        else:
            self.files["contigs"] = os.path.abspath(contigs)
            logger.info(f"Contigs file set as {contigs}")
            self.directory = os.path.dirname(self.files["contigs"])
            logger.info(f"Directory set as {self.directory}")
            self.name = os.path.splitext(os.path.basename(contigs))[0]
            logger.info(f"Name set as {self.name}")

    def load_seqstats(self):
        if not self.files["contigs"]:
            raise Exception(
                "Your Genome doesn't have an input file! Please provide one."
            )

        self.seqstats = dict()

        try:
            stats = subprocess.check_output(
                [CONFIG["third_party"]["seqstats"], self.files["contigs"]]
            ).decode("utf-8")
            stats = stats.split("\n")
            for n in stats:
                n = n.split(":")
                if len(n) == 2:
                    try:
                        self.seqstats[n[0]] = float(n[1].strip().replace(" bp", ""))
                    except (ValueError, IndexError) as error:
                        logger.error(
                            "Something is wrong with the seqstats output. Please check your seqstats command.",
                            exc_info=True,
                        )
            self.seqstats
        except FileNotFoundError:
            raise

    def load_prodigal(
        self, prodigal_out=None, load_geneset=True, load_protset=True, print_=False
    ):
        """
        Attach Prodigal results to class object.

        Input:
        prodigal_out is the output folder generated by the Prodigal function.
        """
        if prodigal_out:
            input = os.path.abspath(prodigal_out)
        else:
            input = self.directory
        # This dictionary is the same as output_files from the prodigal.py script.
        self.files["prodigal"] = dict()

        try:
            os.path.isdir(input)
        except FileNotFoundError:
            print(f"{prodigal_out} directory not found.")

        for file_ in os.listdir(input):
            if file_.startswith(self.name):
                file_ = os.path.join(input, file_)
                if file_.endswith("_genes.fna"):
                    self.files["prodigal"]["genes"] = file_
                elif file_.endswith("_proteins.faa"):
                    self.files["prodigal"]["proteins"] = file_
                elif file_.endswith("_cds.gbk"):
                    self.files["prodigal"]["cds"] = file_
                elif file_.endswith("_scores.txt"):
                    self.files["prodigal"]["scores"] = file_
                else:
                    if print_:
                        logger.debug(
                            f"{file_} apparently is not a Prodigal output file. Ignoring it."
                        )
                    pass
        logger.info(f"Loaded Prodigal files for {self.name}.")

        if load_geneset:
            self.load_geneset()
        if load_protset:
            self.load_protset()

    def load_geneset(self, kind="prodigal", records="dict"):
        """
        Loads gene sets unto Genome.geneset.
        Uses the 'genes' key from the files[kind] dictionary.
        """
        if kind == "prodigal":
            # Origin is the file from which the set came from
            origin = self.files["prodigal"]["genes"]
            try:
                self.geneset["prodigal"] = dict()
                self.geneset["prodigal"]["records"] = get_records(origin, kind=records)
                self.geneset["prodigal"]["origin"] = origin

            except Exception:
                raise
        elif kind in CONFIG["db"].keys():
            origin = self.files[kind]["annotation"]
            try:
                self.geneset[kind] = dict()
                self.geneset[kind]["records"] = get_records(origin, kind=records)
                self.geneset[kind]["origin"] = origin
            except Exception:
                raise
        else:
            logger.error(
                f"Passed {kind} kind of geneset input. Please specify a valid input from either an annotation or Prodigal file.",
                exc_info=True,
            )

        # Maybe change this if/else block later.
        if self.geneset[kind]:
            if records != "gen":
                logger.info(
                    f"Loaded gene set from {kind.capitalize()} data. It has {len(self.geneset[kind]['records'])} genes."
                )
            else:
                logger.info(f"Loaded gene set from {kind.capitalize()} data.")
        else:
            print(f"No gene set found in {origin}.")

    def load_protset(self, kind="prodigal", records="dict"):
        """
        Loads protein sets unto Genome.protset.
        Uses the 'protein' key from the files[kind] dictionary.
        """
        if kind == "prodigal":
            origin = self.files["prodigal"]["proteins"]
            try:
                self.protset["prodigal"] = dict()
                self.protset["prodigal"]["records"] = get_records(origin, kind=records)
                # TODO: attach origin to a variable (stated by 'kind')
                self.protset["prodigal"]["origin"] = origin

            except Exception:
                raise
        elif kind == "prokka":
            origin = self.files["prokka"]["proteins"]
            try:
                self.protset["prokka"] = dict()
                self.protset["prokka"]["records"] = get_records(origin, kind=records)
                self.protset["prokka"]["origin"] = origin
            except Exception:
                raise

        if self.protset[kind]:
            if records == "list":
                logger.info(
                    f"Loaded protein set from {kind.capitalize()} data. It has {len(self.protset[kind]['records'])} proteins."
                )
            else:
                logger.info(f"Loaded protein set from {kind.capitalize()} data.")
        else:
            print(f"No proteins set found in {origin}.")

    """
    Run methods. Possess the @timer_wrapper decorator to measure runtime.
    """

    @timer_wrapper
    def df_prodigal(self, kind="gene"):
        if not self.geneset["prodigal"]:
            try:
                self.load_geneset()
            except Exception("Please load your Prodigal geneset."):
                raise

        self.geneset["prodigal"]["df"] = pd.DataFrame(
            columns=(
                "id",
                "start",
                "stop",
                "strand",
                "id_",
                "partial",
                "start_type",
                "rbs_motif",
                "rbs_spacer",
                "gc_cont",
            )
        )

        def extract_row(record):
            row = record.description.split(";")
            row = row[0].split(" # ") + row[1:]
            row = [i.split("=")[-1] for i in row]
            return row

        for record in self.geneset["prodigal"]["records"]:
            # This one liner appends each record to the end of the dataframe.
            self.geneset["prodigal"]["df"].loc[
                len(self.geneset["prodigal"]["df"])
            ] = extract_row(record)

    @timer_wrapper
    def run_prodigal(self, quiet=True, load_sets=["gene", "prot"]):
        """
        Check for contigs file, run Prodigal on file.
        """
        self.valid_contigs(quiet)
        input = self.files["contigs"]
        logger.info(
            f"Starting Prodigal. Your input file is {input}. Quiet setting is {quiet}."
        )
        p = Prodigal(input, output=self.directory, quiet=quiet)
        self.files["prodigal"] = p.run()
        if "gene" in load_sets:
            self.load_geneset()
        if "prot" in load_sets:
            self.load_protset()

    @timer_wrapper
    def run_prokka(self, quiet=True, load_sets=["gene", "prot"], **kwargs):
        """
        Check for contigs file, run Prokka on file.
        """
        self.valid_contigs(quiet)
        input = self.files["contigs"]
        logger.info(
            f"Starting Prokka. Your input file is {input}. Quiet setting is {quiet}."
        )
        prokka_out = prokka(input, **kwargs)
        self.files["prokka"] = prokka_out
        if "gene" in load_sets:
            self.load_geneset("prokka")
        if "prot" in load_sets:
            self.load_protset("prokka")

    @timer_wrapper
    def blast_seqs(self, db, blast="n", evalue=CONFIG["blast"]["evalue"]):
        """
        Blasts geneset.
        :param db: From config.py db
        :param evalue: evalue to use in Blast
        """
        try:
            db_path = CONFIG["db"][db]
        except KeyError:
            logger.error(
                f"Choose a valid database from {CONFIG['db'][db]}.", exc_info=True
            )
        query = self.files["prodigal"]["genes"]
        out = os.path.join(self.directory, self.name + f"_{db}_blast.xml")
        self.files[db] = dict()
        self.files[db]["xml"] = out
        logger.info(f"Blasting {self.name} to {out}.")

        def blast_method(*args, **kwargs):
            if blast in ("n", "nucl", "nucleotide", "blastn"):
                blast_cmd = NcbiblastnCommandline(*args, **kwargs)
            elif blast in ("p", "prot", "protein", "blastp"):
                blast_cmd = NcbiblastpCommandline(*args, **kwargs)
            elif blast in ("x", "blastx"):
                blast_cmd = NcbiblastxCommandline(*args, **kwargs)
            else:
                raise Exception(
                    "Choose a valid option from 'blastn', 'blastp' or 'blastx'."
                )

            return blast_cmd

        blast_cmd = blast_method(
            query=query,
            db=db_path,
            evalue=evalue,
            out=out,
            outfmt=5,
            num_alignments=5,
            num_threads=CONFIG["threads"],
        )
        stdout, stderr = blast_cmd()
        self.parse_xml_blast(db)

    def parse_xml_blast(self, db, write_hits=True):
        """
        Blasts XML out
        :return:
        """
        hits = []
        with open(self.files[db]["xml"]) as f:
            blast_records = NCBIXML.parse(f)
            for i in blast_records:
                if i.alignments:
                    i.id = i.query.split(" #")[0]
                    i.query = " ".join((i.id, i.alignments[0].hit_def))
                    hits.append(i)
        logger.info(f"Found {len(hits)} hits.\n")

        self.geneset[db] = dict()
        self.geneset[db][
            "origin"
        ] = f"Blast of {self.files['prodigal']['genes']} to {CONFIG['db'][db]}."
        self.geneset[db]["records"] = list()

        for i in hits:
            annotation = self.geneset["prodigal"]["records"][i.id]
            annotation.description = i.query
            self.geneset[db]["records"].append(annotation)
            print(i.query)

        if write_hits:
            out_f = os.path.join(self.directory, self.name + f"_{db}.fasta")
            out_h = os.path.join(self.directory, self.name + f"_{db}.hits")
            self.files[db]["annotation"] = out_f
            self.files[db]["hits"] = out_h
            with open(out_f, "w") as f:
                SeqIO.write(self.geneset[db]["records"], f, "fasta")

            with open(out_h, "w") as f:
                for i in [j.description for j in self.geneset[db]["records"]]:
                    f.write(i + "\n")

            logger.info(
                f"Wrote {len(self.geneset[db]['records'])} annotated sequences to {out_f}."
            )

    def to_json(self, out_path=None):
        """
        Writes a json output of our genome object, that can be load back into Abacat.
        We filter the geneset and protset keys because they would take too much space,
        loading them back when we import the json obj.
        """
        json_out = dict()
        for key, value in self.__dict__.items():
            if key not in ("geneset", "protset"):
                json_out[key] = value

        if not out_path:
            out_path = os.path.join(self.directory, self.name + ".json")
        with open(out_path, "w") as f:
            json.dump(json_out, f, indent=3)
        logger.info(f"Wrote json file of {self.name} to {out_path}.")

    def run_pathways(self, info=True, evalue=10 ** -3):
        """
        Takes the phenotyping geneset records and checks them against the pathways object
        from the CONFIG module.
        :return: pathway genes, a dict containing pathways as keys and identified records as values.
        """
        if "phenotyping" not in self.files.keys():
            logger.info("Phenotyping files not found. Running BLAST now.")
            self.blast_seqs(db="phenotyping", blast="blastx", evalue=evalue)
        elif "phenotyping" not in self.geneset.keys():
            logger.info("Phenotyping records not found. Loading from BLAST out.")
            self.load_geneset(kind="phenotyping")

        self.pathways = dict()
        for (
            k,
            v,
        ) in (
            pathways.items()
        ):  # Note that this is the 'pathways' imported from config.py
            self.pathways[k] = []
            for gene in self.geneset["phenotyping"]["records"]:
                desc = gene.description.split()[1].split(".")[1]
                if desc in v:
                    self.pathways[k].append(gene.description)
            if info:
                logger.info(f"Found {len(self.pathways[k])} genes for {k}.")


@is_fasta_wrapper
def from_fasta(fasta_file, run_prodigal=False, load_prodigal=False):
    """
    A function to load assemblies and run Prodigal directly.

    Input:
    A valid contigs file.

    Returns:
    An Genome object.
    """
    genome = Genome()

    logger.info(f"Loading contigs file from {fasta_file}.")
    genome.load_contigs(fasta_file)
    if run_prodigal:
        print(f"Running Prodigal for {genome.name}.")
        genome.run_prodigal()

    if load_prodigal:
        genome.load_prodigal()

    return genome


def from_json(json_file):
    """
    :param json_file: A json file like the one exported from Genome.to_json()
    :return: A genome instance.
    """
    logger.info(f"Loading genome from {json_file}.")
    with open(json_file) as f:
        j = json.load(f)

    g = Genome(name=j["name"], directory=j["directory"])

    for attribute in ("files", "seqstats", "pathways"):
        if attribute in j.keys():
            setattr(g, attribute, j[attribute])
            logger.info(f"Loaded {attribute} to {g.name}.")

    g.load_geneset()
    g.load_protset()

    del j
    return g
