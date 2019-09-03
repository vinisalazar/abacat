import os
import json

CONFIG = {
    "db": {
        "COG": "/Users/viniWS/Bio/db/COG/prot2003-2014.fa",
        "megares": "/Users/viniWS/Bio/db/megares_v1.01/megares_database_v1.01.fasta",
        "phenotyping": "/Users/viniWS/Bio/db/phenotyping/phenotyping.fasta",
        "pathways": "/Users/viniWS/Bio/bactools/pipelines/pathways.json"
    },
    "threads": int(os.cpu_count() / 2),
    "blast": {"evalue": 10 ** -20},
}


def load_pathways(pathways_json):
    with open(pathways_json) as f:
        pathways = json.load(f)

    return pathways


pathways = load_pathways(CONFIG["db"]["pathways"])
