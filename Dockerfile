FROM python:3.10
RUN apt-get update \
    && apt-get install --no-install-recommends -y libgl1

ENV GIT_SSL_NO_VERIFY=1

# Configure Poetry
ENV POETRY_HOME=/opt/poetry
ENV POETRY_VENV=/opt/poetry-venv
ENV POETRY_CACHE_DIR=/opt/.cache

# Install poetry separated from system interpreter
RUN python3 -m venv $POETRY_VENV \
    && $POETRY_VENV/bin/pip install -U pip setuptools \
    && $POETRY_VENV/bin/pip install poetry

# Add `poetry` to PATH
ENV PATH="${PATH}:${POETRY_VENV}/bin"

WORKDIR /home
RUN git clone https://github.com/ULB-Metronu/georges.git
RUN poetry config virtualenvs.in-project true
WORKDIR /home/georges
RUN poetry install

ENV PATH="/home/georges-core/.venv/bin:$PATH"
RUN mkdir work
WORKDIR /home/work
CMD jupyter-lab --ip 0.0.0.0 --no-browser --allow-root --port=8899
