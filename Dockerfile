FROM continuumio/miniconda3
MAINTAINER Cedric Hernalsteens "cedric.hernalsteens@iba-group.com"
RUN git clone https://github.com/chernals/georges.git
RUN conda env create --file georges/environment.yml
RUN [ "/bin/bash", "-c", ". /opt/conda/etc/profile.d/conda.sh && conda activate ipy3 && pip install ./georges" ]
# To build matplotlib.pyplot font cache
RUN [ "/bin/bash", "-c", ". /opt/conda/etc/profile.d/conda.sh && conda activate ipy3 && ipython -c 'import matplotlib.pyplot'" ]
CMD [ "/bin/bash", "-c", ". /opt/conda/etc/profile.d/conda.sh && conda activate ipy3 && jupyter notebook --ip=0.0.0.0 --port=8888 --allow-root" ]