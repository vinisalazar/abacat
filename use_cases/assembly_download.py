"""
- access downloaded GBK file.
- extract assembly ID.
"""

import os
from Bio import Entrez

input_dir = "/Users/viniWS/storage/neorefs/genomes/bacteria/"

def find_gbk_files(input_dir):
    gbk_files = []
    for dirs, subdirs, files in os.walk(input_dir):
        for file in files:
            if os.path.splitext(file)[1] == ".gbk":
                gbk_files.append(os.path.join(dirs, file))

    return gbk_files

gbk = find_gbk_files(input_dir)


def extract_assembly_accession(list_of_gbk_files, write_to=None):
    acc = []
    for gbk in list_of_gbk_files:
        with open(gbk) as f:
            record = SeqIO.parse(f, "genbank")
            acc.append([i for i in next(record).dbxrefs if "Assembly" in i][0].split("Assembly:")[1])

    if write_to:
        with open(write_to, "w") as f:
            for i in acc:
                f.write(i + "\n")

        if os.path.isfile(write_to):
            print(f"Wrote {len(acc)} accessions to {write_to}.")

    return acc

extract_assembly_accession(gbk, "accessions.txt")
