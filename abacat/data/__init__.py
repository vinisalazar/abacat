"""
Directory for data files. Contains the pathways json file and
also files for the test runs.
"""

from os import path

# This variable returns this directory
data_dir = path.dirname(__file__)
genomes_dir = path.join(data_dir, "genomes")
local_db_dir = path.join(data_dir, "db")
