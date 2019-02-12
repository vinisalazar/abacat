import os
import argparse
import subprocess
import pandas as pd
from Bio import SeqIO

"""
Identify HGT genes from the CompareM output. Extract sequences from .fna or .faa files.
"""

def get_input_files(input_dir):
    """
    Return file names from args.input_dir.
    """
    input_files = os.listdir(input_dir)

    # We only want .csv or .tsv files, with full path
    input_files = [os.path.join(input_dir, i) for i in input_files if i[-2:] == 'sv']

    return input_files


def label_genome_row(input_file, print=True):
    """
    CompareM outputs the genome's row as <complete genome>, we want to label it
    with the organism name.

    We use the sed command with -i (in place).
    """
    if not os.path.isfile(input_file):
        return f"{input_file} not found."

    genome_name = os.path.basename(input_file).split('.genomic')[0]

    os.system(f"sed -i -e 's/\<complete genome\>/{genome_name}/g' {input_file}")
    os.system(f"sed -i -e 's/$/\t{genome_name}/' {input_file}")

    if print:
        print(f"Added genome label to {input_file}")

    return None


def concatenate_input_files(input_dir, output_file, preprocessed=True):

    input_files = get_input_files(input_dir)

    if not preprocessed:
        print(f"Adding labels to {len(input_files)} files.")
        for n in input_files:
            label_genome_row(n, print=False)

    with open(input_files[0]) as f:
        header = f.readline()
        header = "\t".join(header.split('\t')[:-1] + ["Organism\n"])

    print(f"Concatenating {len(input_files)} files. This might take a moment.")
    with open(f"{output_file}", "w") as f:
        f.write(header)
        for file in input_files:
            os.system(f"tail -n +2 {file} >> {output_file}")

    size = str(round(os.path.getsize(output_file) / 2 ** 20, 2)) + " MB"

    print(f"Created {output_file}. Size is {size}.")

    return None



def format_df_filter_hgt(input_file, out_file=None, gc=10.0, bp=300, id_hgt=False):
    """
    Reads an input codon or DI usage file, adds organism and diff_gc columns.

    Can also filter putative HGT genes by the -gc (GC content) and -bp
    (Length) cutoff parameters.
    """
    # Reading the dataframe
    df = pd.read_csv(input_file, sep='\t')
    number_of_genes = len(df)
    print(f"You have {number_of_genes} genes.")  # Keeping the user informed

    # Some wrangling
    df.rename({df.columns[-1]: "Organism"}, inplace=True, axis=1)
    df.GC = df.GC.astype(float)
    df["diff_GC"] = df.apply(lambda row: abs(float(row.GC) - float(df[df["Gene Id"] == row.Organism].GC[0])), axis=1)

    # Select HGT genes, genome name and number of genes.
    # Cutoff params is gc difference from the complete genome and length of the gene. We don't want short genes (< 100-200)
    if id_hgt:
        hgt = df[df.diff_GC.apply(lambda n: n >= gc) & df["Length (bp)"].apply(lambda x: x >= bp)]

        # We want the complete genome at the start of the df
        hgt = hgt.append(df.iloc[0,:]).iloc[::-1]
        number_of_genes = len(hgt) -1

        # Save the resulting dataframe
        out_file = f"{input_file}_hgt_{number_of_genes}g_{gc}_{bp}.tsv"

        if number_of_genes:
            df.to_csv(os.path.join(out_dir, out_file), sep='\t', index=False)
            print(f"Found {number_of_genes} putative HGT genes.\nSaved results to {out_file}.")
        else:
            print(f"Didn't find putative HGT genes for {df_file}. Maybe tweak the params.")

    else: # This is in case we want all the genes.
        if not out_file:
            out_file = f"{genome_name}_codon_usage_format.tsv"
        print(f"Saving {number of genes} to {out_file}.")
        df.to_csv(out_file, sep='\t', index=False)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Takes a directory of CompareM codon usage results and identifies\
                                                  putative HGT genes and parses them to a result file.")

    parser.add_argument("-i", "--input", help="Path to directory of CompareM results.")
    parser.add_argument("-o", "--output", help="Path to output directory.", default='putative_hgt_out/')
    parser.add_argument("-hgt", "--hgt", help="ID HGT genes based on GC and Length values.", default=True)
    parser.add_argument("-gc", "--gc_content", help="GC content difference of gene to complete genome to identify putative HGT genes.", default=15.0)
    parser.add_argument("-bp", "--length", help="Minimum length of gene for it to be considered.", default=300)

    args = parser.parse_args()

    in_dir, out_dir = args.input, args.output
    gc, bp = args.gc_content, args.length



# Usage
# python scripts/id_hgt_genes.py -i ~/Bio/mussismilia/from_mussismilia/lgt_codon_results/ -hgt False -o ~/Bio/masters/codon_usage_format_labels
