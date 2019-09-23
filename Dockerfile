FROM continuumio/miniconda3
RUN sudo apt-get update -y; sudo apt-get upgrade -y
RUN conda install pandas pytest biopython prodigal fastani scipy blast -c conda-forge -c bioconda -y
RUN git clone git@github.com:vinisalazar/abacat.git
RUN cd abacat/
RUN pytest