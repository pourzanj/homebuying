FROM --platform=linux/amd64 rocker/rstudio:4.4.2

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    libmagick++-dev \ 
    && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /homebuying

WORKDIR /homebuying
COPY . /homebuying
#RUN Rscript -e "renv::restore()"
#RUN Rscript -e "cmdstanr::install_cmdstan()"