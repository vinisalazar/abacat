"""
- access downloaded GBK file.
- extract assembly ID.
"""

import os
import argparse
from Bio import SeqIO


def find_gbk_files(input_dir):
    gbk_files = []
    for dirs, subdirs, files in os.walk(input_dir):
        for file in files:
            if os.path.splitext(file)[1] == ".gbk":
                gbk = os.path.join(dirs, file)
            if os.path.splitext(file)[1] == ".fasta":
                fasta = True
            else:
                fasta = False
            if fasta == False:
                gbk_files.append(gbk)
            else:
                pass

    return gbk_files


def extract_assembly_accession(list_of_gbk_files, write_to="accessions.txt"):
    acc = []
    not_found = []
    not_found_p = os.path.splitext(write_to)[0] + "_not_found.txt"
    for gbk in list_of_gbk_files:
        with open(gbk) as f:
            gbk_ = os.path.basename(os.path.splitext(gbk)[0])
            record = SeqIO.parse(f, "genbank")
            try:
                acc_id = [i for i in next(record).dbxrefs if "Genome" in i][0].split(
                    "Genome:"
                )[1]
                acc.append((gbk_, acc_id))
            except IndexError:
                not_found.append(gbk_)
                print(gbk_)
                pass

    if write_to:
        with open(write_to, "w") as f:
            for i in acc:
                f.write(":".join(i) + "\n")

        if len(not_found) >= 1:
            with open(not_found_p, "w") as f:
                for i in not_found:
                    f.write(i)

        if os.path.isfile(write_to):
            print(f"Wrote {len(acc)} accessions to {write_to}.")

        if os.path.isfile(not_found_p):
            print(f"Wrote {len(not_found)} not found {not_found_p}.")

    return acc


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Used to extract Genome accessions from a list of .gbk files inside directories."
    )
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    args = parser.parse_args()

    gbk = find_gbk_files(args.input)
    extract_assembly_accession(gbk, write_to=args.output)
