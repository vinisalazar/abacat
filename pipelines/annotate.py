"""
Example annotation pipeline.

1. Load contigs.

2. Run seqstats.

3. Call genes.

4. Annotate.

"""
import argparse
from bactools import Genome, timer_wrapper

@timer_wrapper
def main(input_, db):
    genome = Genome(input_)
    genome.load_seqstats()
    genome.print_seqstats()
    genome.run_prodigal()
    genome.blastn_seqs(db=db)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotation pipeline. Starts with a contig file.")
    parser.add_argument("-i", "--input", help="Input file. Must a valid FASTA contigs file.")
    parser.add_argument("-db", "--database", type=str, help="Database name. Must be in bactools.CONFIG.py db parameter.")
    args = parser.parse_args()
    main(args.input, args.database)

