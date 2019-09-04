import os
from abacat import timer_wrapper, CONFIG
from Bio import Entrez
from Bio import SeqIO

Entrez.email = CONFIG["email"]


class GenomeFetch:
    """
    Class GenomeFetch. Hosts accession numbers that will be queried.
    """

    def __init__(self, records_file=None):
        super(GenomeFetch, self).__init__()
        self.records_file = records_file
        self.accessions = None
        self.out_dir = CONFIG["out_dir"]
        if records_file:
            if os.path.splitext(records_file)[1] in (
                ".fna",
                ".fasta",
            ):  # Checks if it is a FASTA file.
                self.load_accessions(fasta=True)
            else:
                self.load_accessions()

    def load_accessions(self, fasta=False):
        """
        Load accession numbers as queries.
        Fasta: if True, accession numbers should be in a fasta file as sequence
        descriptions. Else, should be in a file, one per line.
        """
        if fasta:
            with open(self.records_file) as f:
                records = SeqIO.parse(f, "fasta")
                desc = ["_".join(i.description.split("_")[:2]) for i in records]
                self.accessions = [
                    Query(j, out_dir=self.out_dir) for j in {i for i in desc}
                ]  # Don't save the index.
        else:
            with open(self.records_file) as f:
                records = f.readlines()
                records = [i.strip() for i in records]
                self.accessions = [Query(i, out_dir=self.out_dir) for i in records]


class Query(object):
    """
    Query to nucletide database. Should be a valid NCBI accession number.

    It is important to note that some of the accession numbers need to be converted.
    """

    def __init__(self, repr, out_dir=None):
        super(Query, self).__init__()
        self.repr = repr
        self.out_dir = os.path.join(out_dir, repr)
        self.out_gb = os.path.join(self.out_dir, self.repr + ".gbk")
        self.out_fasta = os.path.join(self.out_dir, self.repr + ".fasta")

    def __repr__(self):
        return self.repr

    def fetch_query(self, force=False, gb=True, fasta=True):
        query = self.repr
        out_dir = self.out_dir

        print(f"Searching for query {query}.")

        if not force:
            if os.path.isfile(self.out_gb):
                print(
                    f"{self.out_gb} already exists. Use force or set a different path."
                )
                gb = False
            if os.path.isfile(self.out_fasta):
                print(
                    f"{self.out_fasta} already exists. Use force or set a different path."
                )
                fasta = False

        if not os.path.isdir(self.out_dir):
            os.mkdir(self.out_dir)

        if gb:

            @timer_wrapper
            def write_gb():
                net_handle = Entrez.efetch(
                    db="nucleotide", id=query, rettype="gb", retmode="text"
                )
                r = net_handle.read()
                with open(self.out_gb, "w") as out_handle:
                    out_handle.write(r)
                    net_handle.close()
                    print(f"Saved .gbk record for query {query} at {self.out_gb}.")

            write_gb()

        if fasta:

            @timer_wrapper
            def write_fasta():
                net_handle = Entrez.efetch(
                    db="nucleotide", id=query, rettype="fasta", retmode="text"
                )
                r = net_handle.read()
                if len(r) < 10:
                    print("Failed to find .FASTA for this entry.")
                    return
                with open(self.out_fasta, "w") as out_handle:
                    out_handle.write(r)
                    net_handle.close()
                    print(f"Saved .fasta record for query {query} at {self.out_fasta}.")

            write_fasta()
