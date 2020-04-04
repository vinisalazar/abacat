## 🥑 Abacat - A Bacterial genome Comparison and Annotation Toolkit 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3740294.svg)](https://doi.org/10.5281/zenodo.3740294)

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
conda create -n abacat abacat -c bioconda

# Alternatively, you can install with pip
pip install abacat  

# If you install with pip, remember you'll need the third party packages listed below

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
