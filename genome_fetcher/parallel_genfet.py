"""
This script runs batch GenomeFetch in parallel to fetch queries.
"""

import os
import argparse
from multiprocessing.pool import ThreadPool
from abacat import CONFIG, GenomeFetch, Query
from Bio import Entrez
from operator import methodcaller

Entrez.email = CONFIG["email"]


def main(file_with_accessions):
    """
    file_with_accessions: our input with all accession numbers we want to
    fetch from NCBI.
    """
    genfet = GenomeFetch(file_with_accessions)

    return genfet.accessions


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-nf", "--nofasta", dest="fasta", action="store_false")
    parser.add_argument("-ng", "--nogb", dest="gb", action="store_false")
    parser.add_argument("-f", "--force", dest="force", action="store_true")
    parser.add_argument("-t", "--threads", default=int(os.cpu_count() / 2))
    parser.set_defaults(fasta=True, gb=True, force=False)

    args = parser.parse_args()

    acc = main(args.input)

    with ThreadPool(int(args.threads)) as p:
        f = methodcaller("fetch_query", force=args.force, fasta=args.fasta, gb=args.gb)
        p.map(f, acc)
