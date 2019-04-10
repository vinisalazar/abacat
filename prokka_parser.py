"""
A script to parse Prokka outputs.

Input:
Directory output of Prokka.
"""
import os
from Bio import SeqIO
from helper_functions import is_fasta_wrapper, timer_wrapper


@is_fasta_wrapper
def ffn_parser(ffn_file, write=True):
    """
    Scans Prokka output .ffn file and creates files with SSU seqs,
    known proteins and hypothetical proteins.
    """

    records = SeqIO.parse(ffn_file, "fasta")
    ssu, known_prots, hypothetical_prots = (
        ("SSU", []),
        ("known_prots", []),
        ("hypothetical_prots", []),
    )

    # Generate our output file name.
    output = os.path.splitext(os.path.abspath(ffn_file))[0]

    for seq in records:
        if "16S ribosomal RNA" in seq.description:
            ssu[1].append(seq)
        elif "hypothetical" in seq.description:
            hypothetical_prots[1].append(seq)
        else:
            known_prots[1].append(seq)

    for seq_type, sequences in (ssu, known_prots, hypothetical_prots):
        if len(sequences) >= 1:
            print(f"Found {len(sequences)} {seq_type} sequences in {ffn_file}.")
            if write:
                output_ = output + f"_{seq_type}.fna"
                with open(output_, "w") as f:
                    SeqIO.write(sequences, f, "fasta")
                    print(f"Wrote {len(sequences)} {seq_type} sequences to {output_}.")

    return ssu[1], known_prots[1], hypothetical_prots[1]


# file = "/Users/viniWS/Bio/masters/test_data/own_data/CCMR0080_newscaffold_prokka/PROKKA_04092019.ffn"
# @timer_wrapper
# def main():
#     files = [
#         "/Users/viniWS/Bio/masters/test_data/own_data/CCMR0080_newscaffold_prokka/PROKKA_04092019.ffn",
#         "/Users/viniWS/Bio/masters/test_data/own_data/CCMR0082_newscaffold_prokka/PROKKA_04092019.ffn",
#         "/Users/viniWS/Bio/masters/test_data/own_data/CCMR0083_newscaffold_prokka/PROKKA_04092019.ffn",
#         "/Users/viniWS/Bio/masters/test_data/own_data/CCMR0085_newscaffold_prokka/PROKKA_04092019.ffn",
#     ]
#     for file in files:
#         ssu_finder(file)
#
#
# main()
