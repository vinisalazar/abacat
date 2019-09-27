FROM continuumio/miniconda3
RUN apt-get update -y; apt-get upgrade -y
RUN conda update -n base -c defaults conda -y
RUN conda install pandas pytest biopython matplotlib prodigal fastani scipy blast==2.9 -c conda-forge -c bioconda -y
RUN git clone https://github.com/vinisalazar/abacat.git
WORKDIR /abacat
RUN cd abacat/
RUN pip install .
RUN pytest .
