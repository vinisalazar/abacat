"""
This script runs batch GenomeFetch in parallel to fetch queries.
"""

import os
import argparse
import multiprocessing as mp
from bactools import CONFIG, GenomeFetch, Query
from Bio import Entrez

Entrez.email = CONFIG["email"]

def main(file_with_accessions, fasta=False, gb=True, force=False):
    """
    file_with_accessions: our input with all accession numbers we want to
    fetch from NCBI.
    """
    genfet = GenomeFetch(file_with_accessions)

    with mp.Pool() as p:
        p.map(Query.fetch_query, genfet.accessions)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")

    args = parser.parse_args()

    main(args.input)
