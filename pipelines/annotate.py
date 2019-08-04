"""
Example annotation pipeline.

1. Load contigs.

2. Run seqstats.

3. Call genes.

4. Annotate.

"""

import os
import argparse
from subprocess import Popen, PIPE
from bactools import Assembly, CONFIG, timer_wrapper

@timer_wrapper
def main(input_):
    assembly = Assembly(input_)
    assembly.load_seqstats()
    print(assembly.seqstats)
    assembly.run_prodigal()
    assembly.load_geneset()
    with open(assembly.name + ".blast_out", "w") as f:
        ix = 1
        no = len(assembly.geneset["prodigal"]["records"])
        for gene in assembly.geneset["prodigal"]["records"]:
            id_ = gene.id
            seq = gene.seq
            cmd = f"blastx -query <(echo -e \'>{id_}\\n{seq}\')"
            cmd += f" -db {CONFIG['db']['COG']}"
            cmd += " -evalue 0.001"
            cmd += " -max_target_seqs 1"
            cmd += " -outfmt '6 qseqid sseqid qseq'"

            print(f"Blasting {id_}. Number {ix}/{no}.")
            p = Popen(cmd, stdout=f, stderr=PIPE, shell=True, executable="/bin/bash")
            out, err = p.communicate()
            ix += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotation pipeline. Starts with a contig file.")
    parser.add_argument("-i", "--input", help="Input file. Must a valid FASTA contigs file.")
    args = parser.parse_args()
    main(args.input)

