import os
import json
from pathlib import Path

def db(dir_name, file_name):
    """
    :param dir_name: Directory name inside the $db_path variable
    :param file_name: File name inside $dir_name. Usually a fasta file.
    :return: Joined path of db_path, dir_name, file_name
    """
    db = os.path.join(db_path, dir_name, file_name)
    if not os.path.isfile(db):
        raise Exception(f"Database not found at {db}. Please check if it is a valid path.")
    return db


def load_pathways(pathways_json):
    with open(pathways_json) as f:
        pathways = json.load(f)

    return pathways


db_path = os.path.join(Path.home(), "Bio/db")
CONFIG = {
    "db": {
        "COG": db("COG", "prot2003-2014.fa"),
        "megares": db("megares_v1.01", "megares_database_v1.01.fasta"),
        "phenotyping": db("phenotyping", "phenotyping.fasta"),
        "pathways": "pipelines/pathways.json"
    },
    "threads": int(os.cpu_count() / 2),
    "blast": {"evalue": 10 ** -20},
}
pathways = load_pathways(CONFIG["db"]["pathways"])
