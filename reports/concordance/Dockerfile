FROM debian:buster-slim
LABEL maintainer="https://github.com/pdiakumis"

# misc pkgs
RUN apt-get update && apt-get install -y git wget bash bzip2 zip less && \
    rm -r /var/lib/apt/lists/* && rm -r /var/cache/apt/*

# conda env setup
WORKDIR /home
COPY environment.yml .
ENV CONDA_PREFIX "/home/miniconda"
ENV PATH "${CONDA_PREFIX}/bin:${CONDA_PREFIX}/envs/hailr/bin:${PATH}"
RUN wget --no-verbose https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p ${CONDA_PREFIX} && rm miniconda.sh && \
    conda install --yes -c conda-forge mamba && \
    mamba env create --file environment.yml && mamba clean --yes --all --force-pkgs-dirs

# Concordance scripts + Rmd
ENV SCRIPT_PATH "/home/scripts"
ENV PATH "${SCRIPT_PATH}:${PATH}"
WORKDIR ${SCRIPT_PATH}
COPY concordance concordance.Rmd .
ENV TZ "Australia/Melbourne"

ENTRYPOINT [ "concordance" ]
