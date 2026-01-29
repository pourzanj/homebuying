FROM --platform=linux/amd64 rocker/rstudio:4.4.2

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update && apt-get install -y --no-install-recommends \
    cmake \
    curl \
    libcurl4-openssl-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libmagick++-dev \
    libtbb-dev \
    && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /homebuying

WORKDIR /homebuying

# Copy only renv files first for better layer caching
# This layer only rebuilds when dependency files change
COPY renv.lock renv.lock
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json

# Install R dependencies (expensive operation, cached unless renv.lock changes)
RUN Rscript -e "renv::restore()"

# Install CmdStan (expensive operation, cached unless previous layers change)
RUN Rscript -e "cmdstanr::install_cmdstan()"

# Copy the rest of the application code
# This layer rebuilds on any code change, but doesn't affect the expensive steps above
COPY . /homebuying