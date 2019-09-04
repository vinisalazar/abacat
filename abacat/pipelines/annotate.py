"""
Example annotation pipeline.

1. Load contigs.

2. Run seqstats.

3. Call genes.

4. Annotate.

"""
import argparse
from abacat import Genome, timer_wrapper, CONFIG


def annotate(input_, db, blast, evalue, prodigal=True):
    """
    :param input_: Input file. Must a valid FASTA contigs file (post-assembly).
    :param db: Database name. Must be in abacat.CONFIG.py db parameter.
    :param blast: Blast method. Choose from 'blastn', 'blastp' or 'blastx'. Default is 'blastn'
    :return:
    """
    genome = Genome(input_)
    genome.load_seqstats()
    genome.print_seqstats()
    if prodigal:
        genome.run_prodigal()
    else:
        genome.load_prodigal()
    genome.blast_seqs(db=db, blast=blast, evalue=evalue)
    return genome


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Annotation pipeline. Starts with a contig file."
    )
    parser.add_argument(
        "-i",
        "--input",
        help="Input file. Must a valid FASTA contigs file (post-assembly).",
    )
    parser.add_argument(
        "-db",
        "--database",
        type=str,
        help="Database name. Must be in abacat.CONFIG.py db parameter.",
    )
    parser.add_argument(
        "-b",
        "--blast",
        type=str,
        default="n",
        help="Blast method. Choose from 'blastn', 'blastp' or 'blastx'. Default is 'blastn'",
    )
    parser.add_argument(
        "-e",
        "--evalue",
        type=float,
        default=CONFIG["blast"]["evalue"],
        help="E-value for BLAST. Default is the one set in abacat/config.py",
    )
    args = parser.parse_args()

    @timer_wrapper
    def run():
        annotate(args.input, args.database, args.blast, args.evalue)

    run()
