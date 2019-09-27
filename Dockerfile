# Fetch image
FROM continuumio/miniconda3

# Update ubuntu
RUN apt-get update -y; apt-get upgrade -y

# Update conda and install dependencies
RUN conda update -n base -c defaults conda -y
RUN conda install ipython pandas pytest biopython matplotlib prodigal fastani scipy blast==2.9 -c conda-forge -c bioconda -y

# Install main library
RUN git clone https://github.com/vinisalazar/abacat.git
WORKDIR /abacat
RUN cd abacat/
RUN pip install .

# Run tests
RUN pytest .
