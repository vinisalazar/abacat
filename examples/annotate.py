"""
Example annotation pipeline.

1. Load contigs.

2. Run seqstats.

3. Call genes.

4. Annotate.

"""

from bactools import Assembly

cyano = Assembly("/Users/viniWS/Bio/bactools/data/CCMR0085_newscaffold.fasta")

cyano.seqstats()