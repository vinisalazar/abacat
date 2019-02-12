import os
import argparse
import subprocess
import pandas as pd
#from Bio import SeqIO

"""
Parse and concatenate results from CompareM codon and DI usage files.

Input: directory or file of CompareM results.
Output: Good ol' dataframe.

# TODO:

[ ] Add back GC and length args.
[ ] Fix preprocessed and solo args boolean value.
[ ] Extract sequences from .fna or .faa files.

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

    if print:
        print(f"Added genome label to {input_file}")

    return None


def concatenate_input_files(input_files, output_file, preprocessed=False, hgt=False):

    if not preprocessed:
        print(f"Adding labels to {len(input_files)} files.")
        for n in input_files:
            label_genome_row(n, print=False)

    with open(input_files[0]) as f:
        header = f.readline()
        header = "\t".join(header.split("\t")[:-1] + ["Organism", "diff_GC\n"])

    print(f"Concatenating {len(input_files)} files. Please be patient.")
    with open(f"{output_file}", "w") as f:
        f.write(header)
        for file in input_files:
            df = pd.read_csv(file, sep='\t')
            df.GC = df.GC.astype(float)

            # Column with organism name for when we cat later.
            df["Organism"] = df.iloc[0,0]

            # Difference of GC content from whole genome.
            df["diff_GC"] = df.GC.apply(lambda n: round(abs(n - df.iloc[0, 1])), 2)
            if hgt:
                df = filter_hgt(df)

            df.to_csv(f, header=False, index=False, sep="\t")

    size = str(round(os.path.getsize(output_file) / 2 ** 20, 2)) + " MB"

    print(f"Created {output_file}. Size is {size}.")

    return None


def filter_hgt(df, gc=10.0, bp=300):
    """
    Reads an input codon or DI usage file, adds organism and diff_gc columns.

    Can also filter putative HGT genes by the -gc (GC content) and -bp
    (Length) cutoff parameters.
    """
    print("Filtering HGT genes.")
    hgt = df[df.diff_GC.apply(lambda n: n >= gc) & df["Length (bp)"].apply(lambda x: x >= bp)]
    number_of_genes = len(hgt)
    print(f"Found {number_of_genes} gene in {hgt.iloc[0,0]} ")

    return hgt



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse and concatenate results from CompareM codon and DI usage files. Input: directory or file of CompareM results. Output: Good ol' dataframe.")
    parser.add_argument("-i", "--input", help="Path to input file or input directory of CompareM results.")
    parser.add_argument("-o", "--output", help="Path to output file.", default='codon_usage_parser_out.tsv')
    parser.add_argument("-hgt", "--hgt", help="ID HGT genes based on GC and Length values", default=False)
    parser.add_argument("-gc", "--gc_content", help="GC content difference of gene to complete genome to identify putative HGT genes.", default=15.0)
    parser.add_argument("-bp", "--length", help="Minimum length of gene for it to be considered.", default=300)
    parser.add_argument("-p", "--preprocessed", help="If the input files are preprocessed or not (They should have a row with the organism name", default=False)
    parser.add_argument("-s", "--solo", help="Run for a single file instead of batch. Input should be a file instead of a dir.", default=False)

    args = parser.parse_args()

    if args.solo:
        concatenate_input_files([args.input], args.output, preprocessed=args.preprocessed)
    else:
        concatenate_input_files(get_input_files(args.input), args.output, preprocessed=args.preprocessed, hgt=args.hgt)





# Usage
# python scripts/id_hgt_genes.py -i ~/Bio/mussismilia/from_mussismilia/lgt_codon_results/ -hgt False -o ~/Bio/masters/codon_usage_format_labels
