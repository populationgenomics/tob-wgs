FROM debian:buster-slim
LABEL maintainer="https://github.com/pdiakumis"

# misc pkgs
RUN apt-get update && apt-get install -y git wget bash bzip2 zip less && \
    rm -r /var/lib/apt/lists/* && rm -r /var/cache/apt/*

# conda env setup
WORKDIR /home
ARG CONDA_ENV="hailr"
ARG CONDA_ENV_FILE="conda_env.yml"
COPY $CONDA_ENV_FILE .
ENV CONDA_PREFIX "/home/miniconda"
ENV PATH "${CONDA_PREFIX}/bin:${CONDA_PREFIX}/envs/${CONDA_ENV}/bin:${PATH}"
RUN wget --no-verbose https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p ${CONDA_PREFIX} && rm miniconda.sh && \
    conda install --yes -c conda-forge mamba && \
    mamba env create --file ${CONDA_ENV_FILE} && mamba clean --yes --all --force-pkgs-dirs

# gcs_connector
RUN curl -fSLO https://raw.githubusercontent.com/broadinstitute/install-gcs-connector/master/install_gcs_connector.py && \
    mkdir -p /gsa-key/ && touch /gsa-key/key.json && \
    /home/miniconda/envs/${CONDA_ENV}/bin/python install_gcs_connector.py --key-file /gsa-key/key.json && \
    rm -rf /gsa-key
