"""
Put all assembly report information into a single table.
Usage:

python assembly_report_extractor.py <PATH TO REPORTS> <PATH OF OUTPUT DATAFRAME>

Where <PATH TO REPORTS> is a directory containing assembly reports, and
And <PATH OF OUTPUT DATAFRAME> is your output file name (.csv extension).

# TODO: Fix append method to df.

'Premature optimization is the root of all evil.' Donald Knuth.

"""

import os
import sys
import argparse
import pandas as pd
from abacat.abacat_helper import timer_wrapper


def dict_from_report(report):
    """
    Makes a dictionary from the report, which will then be converted into a dataframe.
    """
    with open(report) as f:
        r = f.readlines()
        r = r[:15]
        r = [line.strip() for line in r]
        r = [line.replace("# ", "") for line in r]
        r = [line.split(":") for line in r]
        r = [(line[0], line[1].lstrip()) for line in r if len(line) > 1]
        dict_ = dict(line for line in r)
        dict_["fname"] = os.path.basename(report).split("_assembly_report.txt")[0]

        return dict_


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A script to put NCBI Genome reports into a csv file."
    )
    parser.add_argument(
        "-i", "--input", help="Genome report or path with assembly_reports."
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path to output .csv file.",
        default="assembly_reports.csv",
    )
    args = parser.parse_args()

    if not args.input:
        parser.print_help()
        sys.exit(0)

    @timer_wrapper
    def main():
        if os.path.isfile(args.input):
            output = os.path.splitext(args.input)[0] + ".csv"
            dict_ = dict_from_report(args.input)
            df = pd.DataFrame(dict_, index=[0])
            if df.shape[1] < 2:
                raise Exception("Invalid input. Please check your input file.")

            df.to_csv(args.output, sep="\t", index=False)
            print(f"Converted {args.input} into {output}")

        elif os.path.isdir(args.input):
            ix = 0
            reports = [
                i for i in os.listdir(args.input) if i.endswith("assembly_report.txt")
            ]
            print(f"You have {len(reports)} reports.")
            for report in reports:
                try:
                    report = os.path.join(os.path.abspath(args.input), report)
                    dict_ = dict_from_report(report)
                    df = pd.DataFrame(dict_, index=[0])
                    if ix == 0:
                        df.to_csv(args.output, sep="\t", index=False)
                    else:
                        with open(args.output, "a") as f:
                            df.to_csv(f, sep="\t", header=False, index=False)
                    ix += 1
                except:
                    pass

            print(f"Done. Reports written to {args.output}.")

    main()
