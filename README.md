## 🥑 Abacat - A Bacterial genome Comparison and Annotation Toolkit 

Abacat (pronounced *"ABBA-cat"*) is a toolkit for working with bacterial whole genome sequencing
(WGS) data. It provides Python objects to represent elements which are common to WGS analysis workflows,
such as coding sequence (CDS) files, containing genes or proteins, alignment methods, or
statistics about your sequences.

It can be used to annotate freshly generated WGS data and parse existing, curated NCBI data. It can also index 
different BLAST or HMM databases for easy alignment searches and produce comparisons using Average Nucleotide Identity 
(ANI) measurements. Lastly, Abacat's classes provides means of tracking *data provenance* in large scale genomics 
experiments, in a way users can monitor and tailor their workflow accordingly.

**Disclaimer: Abacat is in active development and therefore we do not provide any sort of guarantee of its results**.

Contributions are welcome and encouraged!

For an example on how to use Abacat, please look at the [Jupyter tutorial](./tutorial.ipynb). (Note that you will need Jupyter notebook installed.)

## Installing

The easiest way of installing Abacat is using conda:
```
# Create a new env
conda create -n abacat python=3.6

# Install dependencies
conda install biopython prodigal blast fastani -c bioconda

# Clone the package
git clone https://github.com/vinisalazar/abacat.git

# Install the package
cd abacat/
pip install .

# Test the installation
pytest
```


Like most bioinformatics workflows, Abacat has quite a few dependencies. Here are the main ones:

* [Prodigal](https://github.com/hyattpd/Prodigal)
* [FastANI](https://github.com/ParBLiSS/FastANI)
* [Biopython](https://github.com/biopython/biopython)
* [Pandas](https://pandas.pydata.org/)
* [Conda](https://docs.conda.io/en/latest/)
* [Python >= 3.6](https://www.python.org/downloads/)

Consider citing them if you end up using this workflow.
