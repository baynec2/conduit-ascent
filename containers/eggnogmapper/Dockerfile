FROM mambaorg/micromamba:1.5.8

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
RUN micromamba install -y -n base -c conda-forge -c bioconda eggnog-mapper=2.1.12 \
    && micromamba clean --all --yes

# Set work directory
WORKDIR /work

# Define volumes for persistent data
VOLUME ["/opt/eggnog_db", "/work"]

# Default shell
CMD ["/bin/bash"]
