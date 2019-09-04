#### ABACAT - A Bacterial Annotation and Curation Toolkit

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
conda create -n abacat python=3.6

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
`mkdir data/`
`wget https://bit.ly/2ZuZDBa -O data/example_data.fna`

MRSA data:
`wget \
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/021/895/GCF_001021895.1_ASM102189v1/GCF_001021895.1_ASM102189v1_genomic.fna.gz\
 -O data/GCF_001021895.1_ASM102189v1_genomic.fna.gz`


This file contains a scaffolds for a bacterial genome. We can annotate it with:  
`python examples/annotate.py -i data/example_data.fna` 
