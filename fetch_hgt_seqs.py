import os
import argparse
import pandas as pd
from Bio import SeqIO


def process_genes(zip_, seqs_dir, out_dir, blast_type, label=True, verbose=False):

    """
    Extracts reads from genes files.
    zip_ must be a list of tuples containing i[0] = gene and i[1] = seq file name.
    """

    if blast_type.lower() in ("p", "prot", "nr"):
        blast_type = "p"

    genomes = list(set((i[1] for i in zip_)))

    print(f"We have {len(zip_)} genes in {len(genomes)} genomes.")
    print(f"Input directory is {seqs_dir}.")
    print(f"Output directory is {out_dir}.")

    def parse_seq_file(genome):
        if blast_type == "p":
            fname = genome + '.genomic.fna.genes.faa'
        else:
            fname = genome + '.genomic.fna.genes.fna'
        records = SeqIO.parse(os.path.join(seqs_dir, fname), format='fasta')
        return records

    for genome in genomes:
        records = parse_seq_file(genome)
        genes = [i[0] for i in zip_ if i[1] == genome]
        parsed_records = [i for i in records if i.id in genes]
        if label:  # Label the read with the genome name
            for n in parsed_records:
                n.id = n.id + "~" + genome

        if blast_type == "p":
            fname = genome + '_picked_protein.faa'
        else:
            fname = genome + '_picked_genes.fna'

        with open(os.path.join(out_dir, fname), 'w') as f:
            SeqIO.write(parsed_records, f, format="fasta")

        if verbose:
            print(f"Wrote {len(genes)} genes to {os.path.join(out_dir, fname)}.")


# A little helper function for our args. Source: https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Fetch HGT sequences based on\
                                     their contig ID. Input is a dataframe,\
                                     output is directory with seqs.")

    parser.add_argument("-i", "--input_file", help="Dataframe file with gene IDs and organism names.")
    parser.add_argument("-o", "--output_dir", help="Output folder with sequences and stats.")
    parser.add_argument("-s", "--seqs_dir", help="Sequences input directory.")
    parser.add_argument("-b", "--blast_type", help="Blast type. Set 'p' for protein.", default="n")
    parser.add_argument("-gc", "--gc_filter", help="GC content difference filter to apply. Leave 0 for none.", default=10.0, type=float)
    parser.add_argument("-l", "--length", help="Length filter to apply. Leave 0 for none", default=300, type=int)
    parser.add_argument("-v", "--verbose", type=str2bool, nargs='?', const=False, default='false', help="Activate verbose mode.")

    args = parser.parse_args()

    # Checking if input files exist.
    try:
        os.path.isfile(args.input_file) == True
        os.path.isdir(args.seqs_dir) == True
    except FileNotFoundError:
        raise
        print("Your input files weren't found.")

    # Handling output dir
    if not args.output_dir:
        output_dir = os.path.join(os.getcwd(), "fetch_hgt_out/")
    else:
        output_dir = os.path.join(args.output_dir)  
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)

    # Tell the user what's happening
    print(f"Started processing dataframe. Your output path is {output_dir} ")

    df = pd.read_csv(args.input_file, sep="\t")
    print(f"Your dataframe has {len(df)} genes.")

    print(f"Filtering genes with < {args.length} bp.")
    df = df[df["Length (bp)"] >= args.length]

    print(f"Filtering genes with < {args.gc_filter} GC difference.")
    df = df[df["diff_GC"] >= args.gc_filter]

    print(f"Fetching genes.")
    zip_ = list(zip(df["Gene Id"], df["Organism"]))

    process_genes(zip_, args.seqs_dir, output_dir, blast_type=args.blast_type, label=True, verbose=args.verbose)

    print(f"Done. Wrote {len(zip_)} sequences to {output_dir}.")