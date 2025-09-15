# Dockerfile for conduitR
FROM rocker/r-ver:4.4.3

LABEL maintainer="Charlie Bayne <baynec2@gmail.com>"       version="0.0.0.9000"       description="R tools for metaproteomics analysis using conduit and conduit-GUI"

# Install system dependencies
RUN apt-get update -y && apt-get install -y --no-install-recommends         libcurl4-openssl-dev         libxml2-dev         libssl-dev         libfontconfig1-dev         libharfbuzz-dev         libfribidi-dev         fonts-dejavu         libbz2-dev         liblzma-dev         libzstd-dev         libglpk40         git         libgit2-dev         curl &&     apt-get clean &&     rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" &&     R -e "BiocManager::install(c('QFeatures', 'SummarizedExperiment', 'limma', 'Biostrings', 'KEGGREST', 'impute', 'pcaMethods', 'sechm', 'PCAtools'))" &&     R -e "install.packages('remotes', repos='https://cloud.r-project.org')" &&     R -e "remotes::install_github('baynec2/conduitR')"

# Set environment variable
ENV R_LIBS_USER=/usr/local/lib/R/site-library
