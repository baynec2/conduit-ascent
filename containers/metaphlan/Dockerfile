# Base pinned by digest, not the 1.5.8 tag: a tag can be re-pushed, and this
# image's whole point is that a rebuild reproduces the same tools.
FROM mambaorg/micromamba@sha256:475730daef12ff9c0733e70092aeeefdf4c373a584c952dac3f7bdb739601990

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
# Versions pinned to what the in-use image actually contains, verified by
# inspecting /opt/conda/conda-meta in it — not chosen. An unpinned
# `micromamba install` resolves afresh on every build, so two builds of an
# identical Dockerfile produce different tools (see #68); for these packages a
# version change moves results, not just the environment.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
      metaphlan=4.2.4 \
      bowtie2=2.5.5 \
      samtools=1.23.1 \
      pigz=2.8 \
      python=3.13.13 \
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


