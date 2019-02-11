import os
import argparse
import pandas as pd
from Bio import SeqIO

"""
Identify HGT genes from the CompareM output. Extract sequences from .fna or .faa files.
"""

def putative_hgt(df_file, out_dir, gc=15.0, bp=300, id_hgt=True):
    """
    Reads a dataframe, selects putative HGT genes and outputs results to file.
    """
    genome_name = os.path.basename(df_file).split('.genomic')[0]
    bp = int(bp)
    gc = float(gc)

    # Reading the dataframe
    df = pd.read_csv(df_file, sep='\t')
    complete_genome = df.iloc[0,:]

    # Select HGT genes, genome name and number of genes.
    # Cutoff params is gc difference from the complete genome and length of the gene. We don't want short genes (< 100-200)
    if id_hgt == True:
        hgt = df[df.GC.apply(lambda x: abs(x - df.GC.iloc[0]) >= gc) & df["Length (bp)"].apply(lambda x: x >= bp)]
        no_genes = hgt.shape[0]

        # We want the complete genome at the start of the df
        df = hgt.append(complete_genome).iloc[::-1]
        df["diff_GC"] = df.GC.apply(lambda x: round(abs(x - df.GC.iloc[0]), 2))

        # Replace "Complete genome" with the actual genome name
        df.iloc[0,0] = genome_name

        # Save the resulting dataframe
        out_file = f"{genome_name}_hgt_{no_genes}_{int(gc)}_{bp}.tsv"

        if no_genes > 1:
            df.to_csv(os.path.join(out_dir, out_file), sep='\t', index=False)
            print(f"Found {no_genes} putative HGT genes.\nSaved results to {out_file}.")
        else:
            print(f"Didn't find putative HGT genes for {df_file}. Maybe tweak the params.")

    else: # This is in case we want all the genes.
        df.iloc[0,0] = genome_name
        df["Organism"] = genome_name
        df["diff_GC"] = df.GC.apply(lambda x: round(abs(x - df.GC.iloc[0]), 2))
        no_genes = df.shape[0] - 1
        out_file = f"{genome_name}_codon_usage_format.tsv"
        print(f"Found {no_genes} genes. Saving to {os.path.join(out_dir, out_file)}.")
        df.to_csv(os.path.join(out_dir, out_file), sep='\t', index=False)



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

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    no_hgt = []
    for df_file in os.listdir(in_dir):
        df_file = os.path.join(in_dir, df_file)
        try:
            #import pdb; pdb.set_trace()
            putative_hgt(df_file, out_dir, gc=gc, bp=bp)
        except:
            no_hgt.append(df_file)
            print(f"{df_file} is not a valid CompareM LGT results file.")
            pass

    if no_hgt:
        with open(os.path.join(out_dir, 'failed.txt'), 'w') as f:
            for item in no_hgt:
                f.write("%s\n" % item)

# Usage
# python scripts/id_hgt_genes.py -i ~/Bio/mussismilia/from_mussismilia/lgt_codon_results/ -hgt False -o ~/Bio/masters/codon_usage_format_labels
