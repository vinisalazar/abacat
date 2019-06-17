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


def extract_assembly_accession(list_of_gbk_files, write_to="accessions"):
    acc = []
    for gbk in list_of_gbk_files:
        with open(gbk) as f:
            record = SeqIO.parse(f, "genbank")
            acc.append(
                (
                    os.path.splitext(os.path.basename(gbk))[1],
                    [i for i in next(record).dbxrefs if "Assembly" in i][0].split(
                        "Assembly:"
                    )[1],
                )
            )

    if write_to:
        with open(write_to, "w") as f:
            for i in acc:
                f.write(i[1] + "\n")

        with open(write_to + "_paired", "w") as f:
            for i in acc:
                f.write(": ".join(i) + "\n")

        if os.path.isfile(write_to):
            print(f"Wrote {len(acc)} accessions to {write_to}.")

    return acc


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Used to extract Assembly accessions from a list of .gbk files inside directories."
    )
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    args = parser.parse_args()

    gbk = find_gbk_files(args.input)
    extract_assembly_accession(gbk, write_to=args.output)
