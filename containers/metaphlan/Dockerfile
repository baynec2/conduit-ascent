FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="MetaPhlAn Container" \
      org.opencontainers.image.description="Containerized MetaPhlAn (metagenomic profiling) with bowtie2 and samtools via Bioconda" \
      org.opencontainers.image.source="https://bioconda.github.io/recipes/metaphlan/README.html" \
      org.opencontainers.image.licenses="BSD-3-Clause" \
      maintainer="conduit"

SHELL ["/bin/bash", "-lc"]

# Keep root minimal; switch back to non-root for installs
USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates curl bash \
    && rm -rf /var/lib/apt/lists/*

# Install dependencies into base env via micromamba (bioconda + conda-forge)
USER mambauser
ENV MAMBA_DOCKERFILE_ACTIVATE=1
RUN micromamba install -y -n base -c conda-forge -c bioconda \
      metaphlan \
      bowtie2 \
      samtools \
      pigz \
    && micromamba clean -a -y

# Ensure base env binaries are available without needing entrypoint activation
ENV PATH=/opt/conda/bin:${PATH}

# Prepare work and database directories
USER root
RUN mkdir -p /work /db/metaphlan \
    && chown -R mambauser:mambauser /work /db
USER mambauser
WORKDIR /work

# Environment for MetaPhlAn database location
ENV METAPHLAN_DIR=/db/metaphlan
ENV METAPHLAN_BOWTIE2_DB=/db/metaphlan

# Add entrypoint
COPY --chown=mambauser:mambauser entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["metaphlan", "--help"]


