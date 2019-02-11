import os
import argparse
import subprocess


def call_blastn(query, out_dir, outfmt, db):
    print(f"Blasting {query}.")

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    out_file = os.path.join(out_dir, os.path.basename(query.split('genomic.fna')[0])+'blast_out.tsv')

    subprocess.call(
        f"blastn -query {query} -db {db} -outfmt {outfmt} -out {out_file} \
        -remote -num_alignments 20", shell=True
        )

    return f"File created at {out_file}"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Will perform blastn on genes.")
    parser.add_argument("-i", "--in_dir", help="Input directory with genes.")
    parser.add_argument("-o", "--out_dir", help="Output directory for blast outputs.")
    parser.add_argument("-fmt", "--out_fmt", help="Blast format type.", type=int, default=7)
    parser.add_argument("-db", "--ncbi_db", help="NCBI database to search.", default="nucl")

    args = parser.parse_args()

    in_dir = args.in_dir
    files = os.listdir(in_dir)
    files = [os.path.join(in_dir, file) for file in files]

    for query in files:
        try:
            call_blastn(query, args.out_dir, args.out_fmt, args.ncbi_db)
        except:
            print("Error.")
            pass
