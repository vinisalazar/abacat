import abacat
import pytest
import pandas
from os import path

"""
Module for testing the ANIDendrogram class and dependencies.
"""

key_file = path.join(abacat.genomes_dir, "key_file.tsv")
list_of_genomes = list(abacat.CONFIG["test_genomes"].values())
list_of_genomes = [abacat.Genome(file) for file in list_of_genomes]
output_dir = "tmp_test/"
dn = abacat.ANIDendrogram(output_dir=output_dir)

# See if FastANI binary is there
def test_pyani_install():
    assert path.isfile(abacat.CONFIG["third_party"]["fastANI"]), "Please check your FastANI binary path in the Abacat config file. It should match the `which fastANI` command."


def test_ANIDendrogram_class():
    assert type(dn) is abacat.dendrogram.ANIDendrogram


def test_no_genes_files():
    with pytest.raises(Exception):
        assert dn.make_fastani_input(list_of_genomes, kind="instances")


def test_run_prodigal():
    print("Running Prodigal for genomes. This might take a while.")
    ix = 1
    for genome in list_of_genomes:
        print(f"Running for {ix}/{len(list_of_genomes)} genomes.")
        genome.run_prodigal(quiet=True)
        ix += 1


def test_list_of_genomes():
    assert(len(list_of_genomes) == 7)


def test_make_fastani_input():
     dn.make_fastani_input(list_of_genomes=list_of_genomes, kind="instances")
     assert path.isfile(dn.fastani_input)


def test_fastani_run():
    dn.run()
    assert path.isfile(dn.fastani_output)


def test_make_ani_table():
    dn.make_ani_table()
    assert path.isfile(dn.ani_table)


def test_import_key_file():
    dn.import_key_file(key_file)
    assert type(dn.key_df) is pandas.core.frame.DataFrame


def test_make_dendrogram():
    dn.make_dendrogram(color_threshold=5, filter_rename=True)
    assert path.isfile(dn.fig_output + ".png")
