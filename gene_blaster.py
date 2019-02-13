import os
import argparse
import subprocess

"""
Performs BLAST on file or files within a directory.

TODO:
[ ] Add solo/batch mode.
"""


def call_blast(query, out_dir, outfmt, db='nucl', blast_type="n"):
    print(f"Blasting {query}.")

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    out_file = os.path.join(out_dir, os.path.basename(query.split('.')[0])+'_blast_out.tsv')

    if str.lower(blast_type) in ("p", "blastp", "protein"):
        cmd = f"blastp -query {query} -db {db} -outfmt '{outfmt}' -out {out_file} -num_descriptions 1 -num_alignments 5 -remote"
    else:
        cmd = f"blastn -query {query} -db {db} -outfmt {outfmt} -out {out_file} -num_alignments 20"

    subprocess.call(cmd, shell=True)

    return f"File created at {out_file}"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Will perform blast on sequence files inside the input directory.")
    parser.add_argument("-i", "--in_dir", help="Input directory with genes.")
    parser.add_argument("-o", "--out_dir", help="Output directory for blast outputs.")
    parser.add_argument("-fmt", "--out_fmt", help="Blast format type.", default=1)
    parser.add_argument("-db", "--ncbi_db", help="NCBI database to search.")
    parser.add_argument("-b", "--blast_type", help="Blast types, 'n' for blastn and 'p' for blastp.", default="n")

    args = parser.parse_args()

    in_dir = args.in_dir
    files = os.listdir(in_dir)
    files = [os.path.join(in_dir, file) for file in files if file[-3:] in ("fna", "faa", "fasta")]
    print(f"Found {len(files)} files to blast.")

    for query in files:
        try:
            call_blast(query, args.out_dir, args.out_fmt, args.ncbi_db, args.blast_type)
        except:
            raise
            pass
