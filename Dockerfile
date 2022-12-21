FROM python:3.8
RUN apt-get update \
    && apt-get install --no-install-recommends -y libgl1

ENV GIT_SSL_NO_VERIFY=1

# Configure Poetry
ENV POETRY_VERSION=1.2.0
ENV POETRY_HOME=/opt/poetry
ENV POETRY_VENV=/opt/poetry-venv
ENV POETRY_CACHE_DIR=/opt/.cache

# Install poetry separated from system interpreter
RUN python3 -m venv $POETRY_VENV \
    && $POETRY_VENV/bin/pip install -U pip setuptools \
    && $POETRY_VENV/bin/pip install poetry==${POETRY_VERSION}

# Add `poetry` to PATH
ENV PATH="${PATH}:${POETRY_VENV}/bin"

WORKDIR /home
RUN git clone --recurse-submodules --branch develop https://github.com/rtesse/georges.git
RUN poetry config virtualenvs.in-project true
WORKDIR /home/georges-core
RUN poetry install

ENV PATH="/home/georges-core/.venv/bin:$PATH"
RUN mkdir work
WORKDIR /home/work
CMD jupyter-lab --ip 0.0.0.0 --no-browser --allow-root --port=8899
