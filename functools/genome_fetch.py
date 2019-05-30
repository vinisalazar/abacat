import os
from bactools import Assembly, timer_wrapper
from config import CONFIG
from Bio import Entrez
from Bio import SeqIO

Entrez.email = CONFIG["email"]

class GenomeFetch:
    """
    Genome Fetch is the class that will fetch our genomes from NCBI.
    It is attached to a records_file, which contains the terms used for fetching (accession numbers).
    We want to write those files to .gbk and .fasta records.
    The path to those records will be in the CONFIG file.
    """

    def __init__(self, records_file=None):
        super(GenomeFetch, self).__init__()
        self.records_file = records_file
        self.accessions = None
        self.out_gb = CONFIG["gb"]
        self.out_fasta = CONFIG["fasta"]
        if records_file:
            self.load_accessions()

    def load_accessions(self):
        with open(self.records_file) as f:
            records = SeqIO.parse(f, "fasta")
            desc = [i.description for i in records]
        self.accessions = [Query(j, out_gb_dir=self.out_gb, out_fasta_dir=self.out_fasta) for j in {i[:-3] for i in desc}]  # Don't save the index.


class Query(object):
    """docstring for Query."""

    def __init__(self, repr, out_gb_dir=None, out_fasta_dir=None):
        super(Query, self).__init__()
        self.repr = repr
        self.out_gb_dir = None
        self.out_fasta_dir = None
        self.out_gb = os.path.join(out_gb_dir, self.repr + ".gbk")
        self.out_fasta = os.path.join(out_fasta_dir, self.repr + ".fasta")

    def __repr__(self):
        return self.repr

    def fetch_query(self, force=False, gb=True, fasta=True):
        query = self.repr
        out_gb = self.out_gb
        out_fasta = self.out_fasta

        print(f"Searching for query {query}.")

        if not force:
            if os.path.isfile(out_gb):
                print(f"{self.out_gb} already exists. Use force or set a different path.")
                gb = False
            if os.path.isfile(out_fasta):
                print(f"{self.out_fasta} already exists. Use force or set a different path.")
                fasta = False

        if gb:
            @timer_wrapper
            def write_gb():
                net_handle = Entrez.efetch(db="nucleotide", id=query, rettype="gb", retmode="text")
                with open(out_gb, "w") as out_handle:
                    out_handle.write(net_handle.read())
                    net_handle.close()
                    print(f"Saved .gbk record for query {query} at {out_gb}.")

            write_gb()

        if fasta:
            @timer_wrapper
            def write_fasta():
                net_handle = Entrez.efetch(db="nucleotide", id=query, rettype="fasta", retmode="text")
                with open(out_fasta, "w") as out_handle:
                    out_handle.write(net_handle.read())
                    net_handle.close()
                    print(f"Saved .fasta record for query {query} at {out_fasta}.")

            write_fasta()


genfet = GenomeFetch("/Users/viniWS/storage/neorefs/rev6/test/test_rev6.fasta")

assemblies = dict()

for acc in genfet.accessions:
    assemblies[acc.repr] = Assembly(acc.out_fasta)

for k, v in assemblies.items():
    v.run_prodigal()
