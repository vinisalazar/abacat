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
from Bio.Blast.Applications import NcbiblastxCommandline
from bactools import Assembly, CONFIG, timer_wrapper

@timer_wrapper
def main(input_):
    assembly = Assembly(input_)
    assembly.load_seqstats()
    assembly.seqstats
    assembly.run_prodigal()
    # assembly.load_prodigal(assembly.directory)
    assembly.files["prodigal"]["genes"] = "/Users/viniWS/Bio/bactools/data/test_blast_prots.fna"
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

    # blastx = NcbiblastxCommandline(
    #     query=assembly.files["prodigal"]["genes"],
    #     db=CONFIG["db"]["COG"],
    #     out=os.path.join(assembly.directory, "blastx_out"),
    #     num_descriptions=5,
    #     num_alignments=5,
    #     evalue=10**-3
    # )
    # blastx()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotation pipeline. Starts with a contig file.")
    parser.add_argument("-i", "--input", help="Input file. Must a valid FASTA contigs file.")
    args = parser.parse_args()
    main(args.input)

