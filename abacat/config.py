import os
import json
from pathlib import Path

abacat_path = os.path.dirname(__file__)
db_path = os.path.join(Path.home(), "Bio/db")


def db(dir_name, file_name):
    """
    :param dir_name: Directory name inside the $db_path variable
    :param file_name: File name inside $dir_name. Usually a fasta file.
    :return: Joined path of db_path, dir_name, file_name
    """
    db = os.path.join(db_path, dir_name, file_name)
    # if not os.path.isfile(db):
    #    raise Exception(f"Database not found at {db}. Please check if it is a valid path.")
    return db


def load_pathways(pathways_json):
    with open(pathways_json) as f:
        pathways = json.load(f)

    return pathways


CONFIG = {
    "db": {
        "COG": db("COG", "prot2003-2014.fa"),
        "megares": db("megares_v1.01", "megares_database_v1.01.fasta"),
        "phenotyping": db("phenotyping", "phenotyping.fasta"),
        "pathways": os.path.join(abacat_path, "data/pathways.json"),
    },
    "threads": int(os.cpu_count() / 2),
    "blast": {"evalue": 10 ** -20},
    "test_contigs": os.path.join(
        abacat_path, "data/GCF_001021895.1_ASM102189v1_genomic.fna"
    ),
}

pathways = load_pathways(CONFIG["db"]["pathways"])
