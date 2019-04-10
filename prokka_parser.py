"""
A script to parse Prokka outputs.

Our first goal is to extract 16S sequences. (ffn file)

Input:
Directory output of Prokka.
"""
import os
from Bio import SeqIO
from helper_functions import is_fasta_wrapper, timer_wrapper

file = "/Users/viniWS/Bio/masters/test_data/own_data/CCMR0080_newscaffold_prokka/PROKKA_04092019.ffn"


@is_fasta_wrapper
def ssu_finder(input_file, write=True):
    """
    Scans Prokka output .ffn file for 16S sequences.
    """

    records = SeqIO.parse(input_file, "fasta")
    ssu = []
    for seq in records:
        if "16S ribosomal RNA" in seq.description:
            ssu.append(seq)

    if len(ssu) >= 1:
        print(f"Found {len(ssu)} 16S sequences in {input_file}.")
        output = os.path.splitext(os.path.abspath(input_file))[0] + "_ssu.fna"
        if write:
            print(f"Writing {len(ssu)} sequences to {output}.")
            with open(output, "w") as f:
                SeqIO.write(ssu, f, "fasta")
        return ssu
    else:
        print("No SSU sequences found. Please check your file.")


# @timer_wrapper
# def main():
#     files = [
#         "/Users/viniWS/Bio/masters/test_data/own_data/CCMR0080_newscaffold_prokka/PROKKA_04092019.ffn",
#         "/Users/viniWS/Bio/masters/test_data/own_data/CCMR0082_newscaffold_prokka/PROKKA_04092019.ffn",
#         "/Users/viniWS/Bio/masters/test_data/own_data/CCMR0082_newscaffold_prokka/PROKKA_04092019.ffn",
#         "/Users/viniWS/Bio/masters/test_data/own_data/CCMR0082_newscaffold_prokka/PROKKA_04092019.ffn",
#     ]
#     for file in files:
#         ssu_finder(file)
#
#
# main()
