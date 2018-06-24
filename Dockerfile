FROM continuumio/miniconda3
MAINTAINER Cedric Hernalsteens "cedric.hernalsteens@iba-group.com"
RUN git clone https://github.com/chernals/georges.git
RUN conda env create --file environment.yml
RUN conda activate ipy3
# To build matplotlib.pyplot font cache
RUN ipython -c "import matplotlib.pyplot"
CMD start-notebook.sh --NotebookApp.token=''
