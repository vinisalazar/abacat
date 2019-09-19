import abacat
import pytest
from os import path

"""
Module for testing the Genome class and related methods.
"""
input_contigs = abacat.CONFIG["test_contigs"]
g = abacat.Genome()

def test_genome_class():
    """
    :return: asserts if the class instantiates correctly.
    """
    assert type(g) is abacat.genome.Genome


def test_load_contigs():
    """
    :return: asserts if a valid contig loads correctly as a class attribute.
    """
    g.load_contigs(input_contigs)
    errors = []
    assert_values = {
        "name": path.splitext(path.basename(input_contigs))[0],
        "directory": path.dirname(input_contigs),
        "files": {"contigs": path.abspath(input_contigs)}
    }
    for key, value in assert_values.items():
        if getattr(g, key) != value:
            errors.append(key)
    assert not errors, f"Errors in the following attributes: {errors} when invoking Genome.load_contigs()"


def test_invalid_contigs():
    """
    :return: asserts that we return an exception when loading an invalid file as contigs file.
    """
    g = abacat.Genome()
    with pytest.raises(Exception):
        assert g.load_contigs(abacat.CONFIG["db"]["pathways"])


def test_seqstats():
    """
    :return: Tests whether seqstats values return true.
    """
    g.load_contigs(input_contigs)
    g.load_seqstats()
    errors = []
    correct_values = {  # If you change the test_contigs file, modify this.
        'Total n': 2.0,
        'Total seq': 2864278.0,
        'Avg. seq': 1432139.0,
        'Median seq': 1432139.0,
        'N 50': 2839253.0,
        'Min seq': 25025.0,
        'Max seq': 2839253.0
    }
    for (key1, value1), (key2, value2) in zip(g.seqstats.items(), correct_values.items()):
        if key1 != key2 or value1 != value2:
            errors.append(f"{key1} value is {value1} but should be {value2}.")

    assert not errors, f"Errors in the following keys:\n{errors}"


def test_run_prodigal():
    """
    :return: Runs Prodigal for our genome.
    """
    g.load_contigs(input_contigs)
    g.run_prodigal()
    errors = []
    for key, value in g.files["prodigal"].items():
        if not path.isfile(value):
            errors.append(f"{key} file was not found at {value}.")

    assert not errors, f"Errors in the following keys:\n{errors}."


def test_load_prodigal():
    """
    :return: If existing Prodigal data is loaded correctly
    """
    g = abacat.Genome(input_contigs)
    g.load_prodigal()
    errors = []
    for key, value in g.files["prodigal"].items():
        if not path.isfile(value):
            errors.append(f"Couldn't create Prodigal {key} file.")

    assert not errors, f"Errors in the following files:\n{errors}."


def blast_seqs_megares():
    """
    :return: Blasts the Genome object against the Megares database.
    """
    g.blast_seqs("megares")
    errors = []
    for key, value in g.files["megares"].items():
        if not path.isfile(value):
            errors.append(f"Could find {key} file at {value}.")

    assert not errors, f"Errors with the following files:\n{errors}"

