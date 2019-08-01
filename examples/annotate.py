"""
Example annotation pipeline.

1. Load contigs.

2. Run seqstats.

3. Call genes.

4. Annotate.

"""

import os
import argparse
import subprocess
from Bio.Blast.Applications import NcbiblastxCommandline
from bactools import Assembly, CONFIG, timer_wrapper

@timer_wrapper
def main(input_):
    assembly = Assembly(input_)
    assembly.load_seqstats()
    # assembly.run_prodigal()
    assembly.load_prodigal(assembly.directory)
    blastx = NcbiblastxCommandline(
        query=assembly.files["prodigal"]["genes"],
        db=CONFIG["db"]["COG"],
        out=os.path.join(assembly.directory, "blastx_out"),
        num_descriptions=5,
        num_alignments=5,
        evalue=10**-3
    )
    blastx()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotation pipeline. Starts with a contig file.")
    parser.add_argument("-i", "--input", help="Input file. Must a valid FASTA contigs file.")
    args = parser.parse_args()
    main(args.input)

