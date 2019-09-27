# Fetch image
FROM continuumio/miniconda3

# Update ubuntu and install build packages
RUN apt-get update -y; apt-get upgrade -y; apt-get install build-essential -y; apt-get install gcc -y

# Update conda and install dependencies
RUN conda update -n base -c defaults conda -y
RUN conda install ipython pandas pytest biopython matplotlib prodigal fastani scipy blast==2.9 -c conda-forge -c bioconda -y

# Install seqstats from source
RUN cd /root
RUN git clone --recursive https://github.com/clwgg/seqstats.git
RUN make seqstats/
RUN cp seqstats/seqstats /bin
RUN source /root/.bashrc

# Install main library
RUN git clone https://github.com/vinisalazar/abacat.git
WORKDIR /abacat
RUN cd abacat/
RUN pip install .

# Run tests
RUN pytest .
