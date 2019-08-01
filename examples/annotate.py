"""
Example annotation pipeline.

1. Load contigs.

2. Run seqstats.

3. Call genes.

4. Annotate.

"""

import argparse
from bactools import Assembly

def main(input_):
    ably = Assembly(input_)
    ably.load_seqstats()
    ably.run_prodigal()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotation pipeline. Starts with a contig file.")
    parser.add_argument("-i", "--input", help="Input file. Must a valid FASTA contigs file.")
    args = parser.parse_args()
    main(args.input)

