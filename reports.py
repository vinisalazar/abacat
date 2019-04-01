"""
Put all assembly report information into a single table.
Usage:

python assembly_report_extractor.py <PATH TO REPORTS> <PATH OF OUTPUT DATAFRAME>

Where <PATH TO REPORTS> is a directory containing assembly reports, and
And <PATH OF OUTPUT DATAFRAME> is your output file name (.csv extension).
"""

import os
import glob
import argparse
import pandas as pd


def dict_from_report(report):
    """
    report = full path to NCBI Assembly DB assembly report file.
    """
    with open(report) as f:
        r = f.readlines()
        r = r[:15]
        r = [line.strip() for line in r]
        r = [line.replace('# ', '') for line in r]
        r = [line.split(':') for line in r]
        r = [(line[0], line[1].lstrip()) for line in r if len(line) > 1]
        dict_ = dict(line for line in r)
        org = dict_["Organism name"].replace(" (cyanobacteria)", "").replace(" ", "_")
        assembly = dict_["Assembly name"].replace(" ", "_")
        dict_["fname"] = org + "_" + assembly

        return dict_


def get_base_name(report):
    return os.path.basename(report[:-20])


def make_dataframe(dict_, out_name):
    if not out_name:
        out_name = "assembly_reports.csv"
    df = pd.DataFrame(dict_).T
    df['AssemblyFullName'] = df.index
    df.reset_index(level=0, inplace=True, drop=True)
    df.to_csv(out_name, sep='\t')
    return print(f'Table successfuly created at {out_name}.')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A script to put NCBI Assembly reports into a Pandas dataframe.")
    parser.add_argument("-i", "--input", help="Full path to a directory containing assembly reports.")
    parser.add_argument("-o", "--output", help="Path to output .csv file.")
    args = parser.parse_args()

    ls = glob.glob(os.path.join(args.input, '*/*assembly_report.txt'))
    dict_ = dict()
    for report in ls:
        try:
            dict_[get_base_name(report)] = dict_from_report(report)
        except:
            pass

    make_dataframe(dict_, args.output)
