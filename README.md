#### BACTools - Bacterial Assembly Curation Tools

This is is a Python package containing simple command line tools
 for working with bacterial assemblies.

Dependencies are:
* [Prodigal](https://github.com/hyattpd/Prodigal)
* [Biopython](https://github.com/biopython/biopython)
* [Pandas](https://pandas.pydata.org/)
* [Conda](https://docs.conda.io/en/latest/)
* [Python >= 3.6](https://www.python.org/downloads/)


To install:
```
# create a new env
conda create -n bactools python=3.6

# install dependencies
conda install pandas -c conda-forge
conda install biopython prodigal -c bioconda

# install the package
pip install .

# or
python setup.py install
```

#### Examples
To download example data, use [this Google Drive link](https://drive.google.com/file/d/1uEAvYApArhC4lUZ_2i2jpoLM6bg77KeG/view?usp=sharing)
or simply:  
`wget https://bit.ly/2ZuZDBa -O data/example_data.fna`

This file contains a scaffolds for a bacterial genome. We can annotate it with:  
`python examples/annotate.py -i data/example_data.fna` 
