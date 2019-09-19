"""
Find phenotyping genes in a genome. Can be used with either a JSON or contigs input.
"""

import argparse
from abacat import Genome, from_json, timer_wrapper


def main(input_, evalue, json=False):
    if json:
        g = from_json(input_)
    else:
        g = Genome(input_)
        g.run_prodigal()

    g.run_pathways(evalue=evalue)
    g.to_json()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get pathway genes from a genome. Will write json output file."
    )
    parser.add_argument(
        "-i",
        "--input",
        help="Input genome. Must be either: a valid contigs file or a json genome file.",
    )
    parser.add_argument(
        "-e",
        "--evalue",
        type=float,
        default=10 ** -3,
        help="E-value for BLAST to Pathways DB",
    )
    parser.add_argument(
        "-j",
        "--json",
        type=bool,
        default=False,
        help="Specifies that you're using an already processed JSON input.",
    )
    args = parser.parse_args()

    @timer_wrapper
    def run():
        main(args.input, args.evalue, args.json)

    run()
