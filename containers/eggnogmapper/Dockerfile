# Base pinned by digest, not the 1.5.8 tag: a tag can be re-pushed, and this
# image's whole point is that a rebuild reproduces the same tools.
FROM mambaorg/micromamba@sha256:475730daef12ff9c0733e70092aeeefdf4c373a584c952dac3f7bdb739601990

LABEL maintainer="Charlie Bayne <you@example.com>"

ENV DEBIAN_FRONTEND=noninteractive
ENV EGGNOG_DATA_DIR=/opt/eggnog_db
ENV PATH=/opt/conda/bin:$PATH

# Become root to create dirs and install
USER root

# Create directories and fix permissions
RUN mkdir -p /opt/eggnog_db /work \
    && chown -R $MAMBA_USER:$MAMBA_USER /opt/eggnog_db /work

# Switch back to the non-root micromamba user
USER $MAMBA_USER

# Install eggNOG-mapper from bioconda/conda-forge
# Versions pinned to what the in-use image actually contains, verified by
# inspecting /opt/conda/conda-meta in it — not chosen. An unpinned
# `micromamba install` resolves afresh on every build, so two builds of an
# identical Dockerfile produce different tools (see #68); for these packages a
# version change moves results, not just the environment.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
      eggnog-mapper=2.1.12 \
      diamond=2.1.24 \
      hmmer=3.4 \
      python=3.11.15 \
    && micromamba clean --all --yes

# Set work directory
WORKDIR /work

# Define volumes for persistent data
VOLUME ["/opt/eggnog_db", "/work"]

# Default shell
CMD ["/bin/bash"]
