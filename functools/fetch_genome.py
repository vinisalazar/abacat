import os
from config import CONFIG
from Bio import Entrez
from Bio import SeqIO

Entrez.email = "vinicius.salazar@neoprospecta.com"

records = SeqIO.parse(CONFIG["db"], "fasta")
genome = next(records)
query = "_".join(genome.description.split("_")[:-1])

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
