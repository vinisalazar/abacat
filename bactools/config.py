import os

CONFIG = {
    "db": {
        "COG": "/Users/viniWS/Bio/db/COG/prot2003-2014.fa",
        "megares": "/Users/viniWS/Bio/db/megares_v1.01/megares_database_v1.01.fasta",
    },
    "threads": int(os.cpu_count() / 2),
    "blast": {"evalue": 10 ** -20},
}
