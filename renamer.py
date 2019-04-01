import os
import sys
import glob
import gzip
import shutil
import argparse

"""
Rename files in Assembly Directory Structure.
"""


def ls_dirs(batch_dir, suffix='GC*/'):
    """
    Globs a directory looking for assembly directories, that usually start with
    'GCA' or 'GCF'. Other suffixes can be determined with suffix arg.
    """
    list_dirs = glob.glob(os.path.join(batch_dir, suffix))
    return list_dirs


def ls_and_decompress(assembly_dir, unzip=True):
    """
    Globs an assembly directory to get the file names.
    Unzip option decompresses any .gz files in the directory.
    """
    files = glob.glob(os.path.join(assembly_dir, '*'))
    if unzip:
        print(f"Decompressing files in {assembly_dir}.")
        for file in files:
            if file.endswith('.gz'):
                with gzip.open(file, 'r') as f_in, open(file[:-3], 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                os.remove(file)

        # Refresh the 'files' var with unzipped names
        files = glob.glob(os.path.join(assembly_dir, '*'))

    return files


def get_assembly_ids(assembly_dir, gb=True):
    """
    Extracts file from the assembly report file.
    Return tuple containing (filename, assembly_name, organism_name, assembly_accession)
    """
    files = ls_and_decompress(assembly_dir, unzip=False)
    for file in files:
        if file.endswith("assembly_report.txt"):
            assembly_report = file

    if not assembly_report:
        return "Assembly report not found, please check your assembly directory.\n"

    fname = []
    with open(assembly_report) as report:
        r = [line for line in report.readlines()]
        for line in r:
            if line.startswith('# Assembly name:'):
                assembly_name = line.strip().split(':')[1].strip()
                assembly_name = assembly_name.replace(" ", "_")
                fname.append(assembly_name)

            if line.startswith('# Organism name:'):
                organism_name = line.strip().split(':')[1].strip()
                organism_name = organism_name.replace(" (cyanobacteria)", "")
                organism_name = organism_name.replace(" ", "_")
                fname.append(organism_name)

            if line.startswith("# GenBank assembly accession:"):
                assembly_accession = line.strip().split(":")[1].strip()

            if line.startswith("# RefSeq assembly accession:"):
                assembly_accession = line.strip().split(":")[1].strip()
                if gb:
                    assembly_accession = assembly_accession.replace('F', 'A')

            elif 'assembly accession' in line:
                assembly_accession = line.strip().split(":")[1].strip()

        # The organism name comes first.
        fname = fname[1] + '_' + fname[0]

        return fname, assembly_name, organism_name, assembly_accession


def rename_assembly_dir(assembly_dir):
    """
    Renames files in assembly dir.
    """
    print(f"Renaming files in {assembly_dir}.\n")
    assembly_dir = os.path.abspath(assembly_dir)
    fname, ass_name, org_name, ass_acc = get_assembly_ids(assembly_dir)
    preffix = '_'.join([ass_acc, ass_name])
    files = ls_and_decompress(assembly_dir, unzip=False)

    for file in files:
        if preffix in os.path.basename(file):
            new_file = os.path.basename(file).replace(preffix, fname)
            os.rename(file, os.path.join(assembly_dir, new_file))
        else:
            new_file = '_'.join([fname, os.path.basename(file)])
            os.rename(file, os.path.join(assembly_dir, new_file))

    os.rename(assembly_dir, os.path.join(os.path.dirname(assembly_dir),
        org_name + "_" + ass_name))

    # TODO: Improve this return statement.
    return "Done."


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Renames assembly directories.")
    parser.add_argument("-i", "--input", help="A single assembly directory to be renamed.")
    parser.add_argument("-b", "--batch", help="Directory containing multiple directories to be renamed.")
    # TODO: Add verbose argument.
    args = parser.parse_args()
    if args.batch:
        batch_dir = args.batch
        ls_dirs = ls_dirs(batch_dir)

        for directory in ls_dirs:
            try:
                ls_and_decompress(directory)
                rename_assembly_dir(directory)
            except:
                print(f"Error in {directory}.\n")

    elif args.input:
        dir = args.input
        ls_and_decompress(dir)
        rename_assembly_dir(dir)

    else:
        parser.print_help()
