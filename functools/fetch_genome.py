import os
from config import CONFIG
from Bio import Entrez
from Bio import SeqIO

Entrez.email = "vinicius.salazar@neoprospecta.com"

records = SeqIO.parse(CONFIG["db"], "fasta")
genome = next(records)
query = "_".join(genome.description.split("_")[:-1])


class Genome_Fetch:
    """
    Genome Fetch is the class that will fetch our genomes from NCBI.
    It is attached to a records_file, which contains the terms used for fetching (accession numbers).
    We want to write those files to .gbk and .fasta records.
    The path to those records will be in the CONFIG file.
    """

    def __init__(self, records_file=None):
        super(Genome_Fetch, self).__init__()
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
        self.accessions = {i[:-3] for i in desc}  # Don't save the index.

    def fetch_query(query, force=False, gb=True, fasta=True):
        out_gb = os.path.join(self.out_gb, query + ".gbk")
        out_fasta = os.path.join(self.out_fasta, query + ".fasta")
        if not force:
            if os.path.isfile(out_gb):
                return f"{out_gb} already exists. Use force or set a different path."
            if os.path.isfile(out_fasta):
                return f"{out_fasta} already exists. Use force or set a different path."

        print(f"Searching for query {query}.")

        if gb:
            net_handle = Entrez.efetch(db="nucleotide", id=query, rettype="gb", retmode="text")
            with open(out_gb, "w") as out_handle:
                out_handle.write(net_handle.read())
                net_handle.close()
                print(f"Saved .gbk record for query {query} at {out_gb}.")

        if fasta:
            net_handle = Entrez.efetch(db="nucleotide", id=query, rettype="fasta", retmode="text")
            with open(out_fasta, "w") as out_handle:
                out_handle.write(net_handle.read())
                net_handle.close()
                print(f"Saved .fasta record for query {query} at {out_fasta}.")






def download_query(query, force=False, fasta=True, gb=True):
    out_gb = os.path.join(CONFIG["gb"], query + ".gbk")
    out_fasta = os.path.join(CONFIG["fasta"], query + ".fasta")
    if not force:
        if os.path.isfile(out_gb):
            return f"{out_gb} already exists. Use force or set a different path."
        if os.path.isfile(out_fasta):
            return f"{out_gb} already exists. Use force or set a different path."

    print(f"Searching for query {query}.")

    if gb:
        net_handle = Entrez.efetch(db="nucleotide", id=query, rettype="gb", retmode="text")
        with open(out_gb, "w") as out_handle:
            out_handle.write(net_handle.read())
            net_handle.close()
            print(f"Saved .gbk record for query {query} at {out_gb}.")

    if fasta:
        net_handle = Entrez.efetch(db="nucleotide", id=query, rettype="fasta", retmode="text")
        with open(out_fasta, "w") as out_handle:
            out_handle.write(net_handle.read())
            net_handle.close()
            print(f"Saved .fasta record for query {query} at {out_fasta}.")

download_query(query, force=True)
