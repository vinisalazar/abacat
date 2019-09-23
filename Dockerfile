FROM continuumio/miniconda3
RUN apt-get update -y; apt-get upgrade -y
RUN conda install pandas pytest biopython prodigal fastani scipy blast==2.9 -c conda-forge -c bioconda -y
RUN git clone https://github.com/vinisalazar/abacat.git
RUN cd abacat/
RUN pytest