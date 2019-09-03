"""
Example annotation pipeline.

1. Load contigs.

2. Run seqstats.

3. Call genes.

4. Annotate.

"""
import argparse
from bactools import Genome, timer_wrapper


def annotate(input_, db, blast):
    """
    :param input_: Input file. Must a valid FASTA contigs file (post-assembly).
    :param db: Database name. Must be in bactools.CONFIG.py db parameter.
    :param blast: Blast method. Choose from 'blastn', 'blastp' or 'blastx'. Default is 'blastn'
    :return:
    """
    genome = Genome(input_)
    genome.load_seqstats()
    genome.print_seqstats()
    genome.run_prodigal()
    genome.blastn_seqs(db=db, blast=blast)
    return genome


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotation pipeline. Starts with a contig file.")
    parser.add_argument("-i", "--input", help="Input file. Must a valid FASTA contigs file (post-assembly).")
    parser.add_argument("-db", "--database", type=str, help="Database name. Must be in bactools.CONFIG.py db parameter.")
    parser.add_argument("-b", "--blast", type=str, default="n", help="Blast method. Choose from 'blastn', 'blastp' or 'blastx'. Default is 'blastn'")
    args = parser.parse_args()

    @timer_wrapper
    def run():
        annotate(args.input, args.database, args.blast)

    run()
