FROM ubuntu:20.04

RUN apt-get update \
    && apt-get install --no-install-recommends -y \
        wget \
        build-essential \
        ca-certificates \
        curl \
        git \
        libgl1 \
        python3 \
        python3-pip \
        python3-distutils

RUN useradd app \
    && mkdir -p /home/app \
    && chown -v -R app:app /home/app

RUN curl -sSL https://install.python-poetry.org | python3 -
ENV PATH="/root/.local/bin:$PATH"

WORKDIR /home/app
ENV GIT_SSL_NO_VERIFY=1

RUN git clone --recurse-submodules --branch develop https://github.com/ULB-Metronu/georges.git
RUN poetry config virtualenvs.in-project true
WORKDIR /home/app/georges
RUN poetry install -E sphinx
ENV PATH="/home/app/georges/.venv/bin:$PATH"

# Run test
RUN pytest tests/

WORKDIR /home/app
RUN mkdir reps
WORKDIR /home/app/reps