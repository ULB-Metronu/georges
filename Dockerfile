FROM jupyter/scipy-notebook
# https://github.com/jupyter/docker-stacks/tree/master/scipy-notebook
MAINTAINER Cedric Hernalsteens "cedric.hernalsteens@iba-group.com"
USER root
RUN wget http://madx.web.cern.ch/madx/releases/last-dev/madx-linux64-gnu -O /usr/local/bin/madx
RUN chmod +x /usr/local/bin/madx
USER $NB_USER
RUN git clone https://github.com/chernals/georges.git packages/georges
ENV PYTHONPATH "/home/$NB_USER/work/packages:$PYTHONPATH"
# To build matplotlib.pyplot font cache
RUN ipython -c "import matplotlib.pyplot"
CMD start-notebook.sh --NotebookApp.token=''
