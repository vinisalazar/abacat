import os
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from matplotlib import pyplot as plt

# This data contains whole genome codon usage followed by the genomes identified genes.
df = pd.read_csv("~/Bio/masters/codon_usage_format_labels/codon_usage_all_header.tsv", sep="\t")

# We don't want very short fragments (>= 300)
df_ = df[df["Length (bp)"] >= 300]

# Let's flag putative HGT genes (different GC content than the whole genome)
hgt = df_[df_.diff_GC >= 15.0]

# Some plots
#sns.distplot(hgt["Manhattan distance"])
#sns.regplot(hgt.diff_GC, hgt["Manhattan distance"])
#sns.jointplot(hgt.diff_GC, hgt["Length (bp)"])

# Let's hunt down the seqs
seqs_dir = "/Users/viniWS/Bio/mussismilia/from_mussismilia/genes_predicted"
out_dir = "/Users/viniWS/Bio/masters/putative_hgt_seqs"

zip_ = list(zip(hgt["Gene Id"], hgt.Organism))


def process_genes(zip_, seqs_dir, out_dir, label=True):

    """
    Extracts reads from genes files.
    zip_ must be a list of tuples containing i[0] = gene and i[1] = seq file name.
    """

    genomes = list(set((i[1] for i in zip_)))
    print(f"We have {len(zip_)} genes in {len(genomes)} genomes.")
    print(f"Input directory is {seqs_dir}.")
    print(f"Output directory is {out_dir}.")

    def parse_fna_file(genome):
        fname = genome + '.genomic.fna.genes.fna'
        records = SeqIO.parse(os.path.join(seqs_dir, fname), format='fasta')
        return records

    for genome in genomes:
        records = parse_fna_file(genome)
        genes = [i[0] for i in zip_ if i[1] == genome]
        parsed_records = [i for i in records if i.id in genes]
        if label: # Label the read with the genome name
            for n in parsed_records:
                n.id = n.id + "-from-" + genome

        fname = genome + '_picked_genes.fna'
        with open(os.path.join(out_dir, fname), 'w') as f:
            SeqIO.write(parsed_records, f, format="fasta")

        print(f"Wrote {len(genes)} genes to {os.path.join(out_dir, fname)}.")


def parse_picked_genes(in_dir, out_dir):
    """
    Concatenate output of process_genes() into a single fasta file.
    Seq record must contain genome name as well.

    This is only needed if label = False in process_genes(). Else, we can use !cat.
    """
    files = os.listdir(in_dir)

    for file in files:
        genome_name = file.split("_picked_genes.fna")[0]
        records = list(SeqIO.parse(os.path.join(in_dir, file), format="fasta"))
        for rec in records:
            rec.id = rec.id + "-from-" + genome_name

        out_file = os.path.join(out_dir, genome_name + "_picked_genes_format.fna")

        with open(out_file, "w") as f:
            SeqIO.write(records, f, format="fasta")

        print(f"Wrote formatted file to {out_file}.")


# in_dir = "/Users/viniWS/Bio/masters/putative_hgt_seqs"
# out_dir = "/Users/viniWS/Bio/masters/putative_hgt_seqs_format"
# parse_picked_genes(in_dir, out_dir)
