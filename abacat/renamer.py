import argparse
import gzip
import os
import shutil
import sys
from abacat.abacat_helper import timer_wrapper

"""
A script to rename files in assembly directory structure.

Input:
One or more assembly directories.
Output:
One or more assembly directories with the directory and all files renamed.
"""


def ls_and_decompress(assembly_dir, unzip=True):
    """
    This function lists and decompresses the files in the assembly directory.
    Must be in NCBI ftp file hierarchy. e.g.
    ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/770/275/GCF_000770275.1_ASM77027v1

    Input:
    Genome directory.
    Returns:
    List of files, decompressed.
    """
    # I prefer dealing with absolute paths.
    assembly_dir = os.path.abspath(assembly_dir)

    # Quick error handling
    if not os.path.isdir(assembly_dir):
        raise FileNotFoundError

    files = [os.path.join(assembly_dir, file) for file in os.listdir(assembly_dir)]

    if unzip:
        print(f"Decompressing files in {assembly_dir}.")
        for file in files:
            if file.endswith(".gz"):
                with gzip.open(file, "r") as f_in, open(file[:-3], "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
                os.remove(file)

        # Refresh the 'files' var with unzipped names
        files = [os.path.join(assembly_dir, file) for file in os.listdir(assembly_dir)]

    return files


def parse_assembly_report(assembly_dir):
    """
    Parses information from assembly_report. Uses information to rename.

    Input:
    Genome directory.

    Returns:
    Tuple containing (organism_name, assembly_name, assembly_accession)

    TODO: Improve this function. See RefSeq/GenBank use cases. Major error handling.
    """
    files = ls_and_decompress(assembly_dir, unzip=False)

    for file in files:
        if file.endswith("assembly_report.txt"):
            assembly_report = file

    if not assembly_report:
        return "Genome report not found, please check your assembly directory.\n"

    with open(assembly_report) as report:
        r = [line for line in report.readlines()]
        for line in r:
            if line.startswith("# Genome name:"):
                assembly_name = line.strip().split(":")[1].strip()

            if line.startswith("# Organism name:"):
                organism_name = "_".join(line.strip().split(":")[1].strip().split())

            if line.startswith("# GenBank assembly accession:"):
                assembly_accession = line.strip().split(":")[1].strip()
                if "GCF" in assembly_report:
                    assembly_accession = assembly_accession.replace("GCA", "GCF")

    return assembly_name, organism_name, assembly_accession


def rename_assembly(assembly_dir, rename="organism assembly", parse_organism_name=3):
    """
    Use information from the assembly report to rename our directory and files.

    Input:
    Genome directory.

    Returns:
    Renamed assembly directory and files.

    rename:
    How the file should be renamed.

    parse_organism_name INT:
    Length of organism name epithets.

    Todo: Add renaming options, like assembly + accession, organism + accession.
          Maybe write a helper function to do this.
    """
    assembly_dir = os.path.abspath(assembly_dir)

    print(f"Renaming files in {assembly_dir}\n")

    assembly_name, organism_name, assembly_accession = parse_assembly_report(
        assembly_dir
    )
    preffix = "_".join([assembly_accession, assembly_name])
    files = ls_and_decompress(assembly_dir, unzip=False)

    # This code is unnecessary for now, but it will help with the renaming options later
    if "organism" in rename:
        # This if defines the length of the organism epithets.
        # For 3 we get Genus, species and strain (default)
        if len(organism_name.split("_")) > parse_organism_name:
            organism_name = organism_name.split("_")[:parse_organism_name]
            organism_name = "_".join(organism_name)

        fname = "_".join([organism_name, assembly_name])

    # These characters will be removed from the final filename.
    # This is because they might interefere with posix path.
    bad_chars = "\\/.?()[]{}#<>º´'"
    for char in bad_chars:
        fname = fname.replace(char, "")

    # Rename the files
    for file in files:
        if preffix in os.path.basename(file):
            new_file = os.path.basename(file).replace(preffix, fname)
            os.rename(file, os.path.join(assembly_dir, new_file))

    # Rename the directory
    os.rename(assembly_dir, os.path.join(os.path.dirname(assembly_dir), fname))
    return os.path.join(os.path.dirname(assembly_dir), fname)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Renames assembly directories.")
    parser.add_argument(
        "-i",
        "--input",
        help="Genome directory or directory containing multiple assembly directories.",
    )
    args = parser.parse_args()

    if not args.input:
        parser.print_help()
        sys.exit(0)

    @timer_wrapper
    def main():
        # Check if it is a single assembly directory or multiple:
        if "annotation_hashes.txt" in os.listdir(args.input):
            ls_and_decompress(args.input)
            rename_assembly(args.input)

        else:
            directories = os.listdir(args.input)
            success, failure = 0, 0
            for directory in directories:
                try:
                    directory = os.path.join(args.input, directory)
                    ls_and_decompress(directory)
                    new_dir = rename_assembly(directory)
                    if os.path.isdir(new_dir):
                        success += 1

                except Exception:
                    print(
                        f"Something went wrong with {directory}. Please confirm if it is a proper assembly directory."
                    )
                    failure += 1
                    pass

            print(
                f"Done. You successfully renamed {success} directories and had {failure} failures."
            )

    main()
