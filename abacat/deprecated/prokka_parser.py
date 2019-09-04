"""
A script to parse Prokka outputs.

Input:
Directory output of Prokka.
"""
import argparse
import os
import sys
from Bio import SeqIO
from abacat.abacat_helper import is_fasta_wrapper, timer_wrapper


@is_fasta_wrapper
def ffn_parser(fasta_file, write=True):
    """
    Scans Prokka output .ffn file and creates files with SSU seqs,
    known proteins and hypothetical proteins.

    # TODO: add RNA support (tRNA, ribonucleases, etc.)
    """

    records = SeqIO.parse(fasta_file, "fasta")
    ssu, known_prots, hypothetical_prots = (
        ("SSU", []),
        ("known_prots", []),
        ("hypothetical_prots", []),
    )

    # Generate our output file name.
    output = os.path.splitext(os.path.abspath(fasta_file))[0]

    for seq in records:
        if "16S ribosomal RNA" in seq.description:
            ssu[1].append(seq)
        elif "hypothetical" in seq.description:
            hypothetical_prots[1].append(seq)
        else:
            known_prots[1].append(seq)

    for seq_type, sequences in (ssu, known_prots, hypothetical_prots):
        if len(sequences) >= 1:
            print(f"Found {len(sequences)} {seq_type} sequences in {fasta_file}.")
            if write:
                output_ = output + f"_{seq_type}.fna"
                with open(output_, "w") as f:
                    SeqIO.write(sequences, f, "fasta")
                    print(f"Wrote {len(sequences)} {seq_type} sequences to {output_}.")

    return ssu[1], known_prots[1], hypothetical_prots[1]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A script to parse Prokka output files."
    )
    parser.add_argument(
        "-i", "--input", help="Directory contaning Prokka output files."
    )

    args = parser.parse_args()

    if not args.input:
        parser.print_help()
        sys.exit(0)

    @timer_wrapper
    def main():
        if not os.path.isdir(args.input):
            raise Exception("Your directory wasn't found.")

        for file in os.listdir(args.input):
            file = os.path.join(os.path.abspath(args.input), file)
            if file.endswith(".ffn"):
                ffn_parser(file)

    main()
