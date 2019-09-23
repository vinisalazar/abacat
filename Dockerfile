FROM continuumio/miniconda3
RUN sudo apt-get upgrade; sudo apt-get update
RUN conda install pandas pytest biopython prodigal fastani scipy blast -c conda-forge -c bioconda
RUN git clone git@github.com:vinisalazar/abacat.git
RUN cd abacat/
RUN pytest