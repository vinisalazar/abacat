import os
import json
import subprocess
from pathlib import Path
from abacat.data import data_dir, genomes_dir, local_db_dir

abacat_path = os.path.dirname(__file__)
external_db_path = os.path.join(Path.home(), "Bio/db")


def db(dir_name, file_name, external=False):
    """
    :param dir_name: Directory name inside the $db_path variable
    :param file_name: File name inside $dir_name. Usually a fasta file.
    :param external: True if the database is not packaged within Abacat
    :return: Joined path of db_path, dir_name, file_name
    """
    if not external:
        db = os.path.join(local_db_dir, dir_name, file_name)
    else:
        db = os.path.join(external_db_path, dir_name, file_name)
    # if not os.path.isfile(db):
    #    raise Exception(f"Database not found at {db}. Please check if it is a valid path.")
    return db


def load_phenotyping(phenotyping_json):
    with open(phenotyping_json) as f:
        phenotyping = json.load(f)

    return phenotyping


def get_third_party_bins(binaries=("fastANI", "seqstats")):
    """
    Find path of third party binaries
    :param binaries: List or tuple with name of binaries.
    :return: dict with keys as the name of the binaries, values as the path in the system.
    """
    bins = dict()
    def get_bin(bin_name):
        out = subprocess.Popen(["which", bin_name], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout, stderr = out.communicate()
        return Path(stdout.strip().decode("utf-8"))

    for bin in binaries:
        bins[bin] = get_bin(bin)

    return bins


CONFIG = {
    "db": {
        "COG": db("COG", "prot2003-2014.fa", external=True),
        "megares": db("megares", "megares_database_v1.01.fasta"),
        "phenotyping": db("phenotyping", "phenotyping.fasta"),
        "pathways": db("phenotyping", "pathways.json"),
    },
    "third_party": get_third_party_bins(),  # Docker config
    "threads": int(os.cpu_count() / 2),
    "blast": {"evalue": 10 ** -20},
    "test_genomes": {
        "Staphylococcus aureus CA15": os.path.join(
            genomes_dir, "GCF_001021895.1_ASM102189v1_genomic.fna"
        ),
        "Staphylococcus aureus KUN1163": os.path.join(
            genomes_dir, "GCF_008619075.1_ASM861907v1_genomic.fna"
        ),
        "Staphylococcus epidermidis ATCC 12228": os.path.join(
            genomes_dir, "GCF_000007645.1_ASM764v1_genomic.fna"
        ),
        "Synechococcus elongatus PCC 6301": os.path.join(
            genomes_dir, "GCF_000010065.1_ASM1006v1_genomic.fna"
        ),
        "Synechococcus elongatus PCC 7942": os.path.join(
            genomes_dir, "GCF_000012525.1_ASM1252v1_genomic.fna"
        ),
        "Prochloroccocus marinus MIT 9313": os.path.join(
            genomes_dir, "GCF_000011485.1_ASM1148v1_genomic.fna"
        ),
        "Escherichia coli ATCC 1175": os.path.join(
            genomes_dir, "GCF_003697165.2_ASM369716v2_genomic.fna"
        ),
    },
    "data_dir": data_dir,
}

pathways = load_phenotyping(CONFIG["db"]["pathways"])
