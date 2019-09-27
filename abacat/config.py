import os
import json
from pathlib import Path
from abacat.data import data_dir

abacat_path = os.path.dirname(__file__)

external_db_path = os.path.join(Path.home(), "Bio/db")
local_db_path = os.path.join(abacat_path, "data/db/")


def db(dir_name, file_name, local=False):
    """
    :param dir_name: Directory name inside the $db_path variable
    :param file_name: File name inside $dir_name. Usually a fasta file.
    :return: Joined path of db_path, dir_name, file_name
    """
    if local:
        db = os.path.join(local_db_path, dir_name, file_name)
    else:
        db = os.path.join(external_db_path, dir_name, file_name)
    # if not os.path.isfile(db):
    #    raise Exception(f"Database not found at {db}. Please check if it is a valid path.")
    return db


def load_phenotyping(phenotyping_json):
    with open(phenotyping_json) as f:
        phenotyping = json.load(f)

    return phenotyping


CONFIG = {
    "db": {
        "COG": db("COG", "prot2003-2014.fa", local=False),
        "megares": db("megares_v1.01", "megares_database_v1.01.fasta", local=True),
        "phenotyping": db("phenotyping", "phenotyping.fasta", local=True),
        "phenotyping": os.path.join(data_dir, "phenotyping.json"),
    },
    "third_party": {"fastANI": "/Users/viniWS/Bio/FastANI/fastANI"},
    "threads": int(os.cpu_count() / 2),
    "blast": {"evalue": 10 ** -20},
    "test_contigs": os.path.join(data_dir, "GCF_001021895.1_ASM102189v1_genomic.fna"),
    "data_dir": data_dir,
}

phenotyping = load_phenotyping(CONFIG["db"]["phenotyping"])
