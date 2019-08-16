"""
Example annotation pipeline.

1. Load contigs.

2. Run seqstats.

3. Call genes.

4. Annotate.

"""
import argparse
from bactools import Assembly, timer_wrapper

@timer_wrapper
def main(input_, db):
    assembly = Assembly(input_)
    assembly.load_seqstats()
    assembly.print_seqstats()
    assembly.run_prodigal()
    assembly.blastn_seqs(db=db)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotation pipeline. Starts with a contig file.")
    parser.add_argument("-i", "--input", help="Input file. Must a valid FASTA contigs file.")
    parser.add_argument("-db", "--database", type=str, help="Database name. Must be in bactools.CONFIG.py db parameter.")
    args = parser.parse_args()
    main(args.input, args.database)

