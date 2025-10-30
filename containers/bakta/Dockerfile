FROM mambaorg/micromamba:1.5.8

LABEL maintainer="Charlie Bayne <you@example.com>"

ENV DEBIAN_FRONTEND=noninteractive
ENV BAKTA_DB=/opt/bakta_db
ENV PATH=/opt/conda/bin:$PATH

# Become root to create dirs and install
USER root

# Create directories and fix permissions
RUN mkdir -p /opt/bakta_db /work \
    && chown -R $MAMBA_USER:$MAMBA_USER /opt/bakta_db /work

# Switch back to the non-root micromamba user
USER $MAMBA_USER

# Install Bakta from bioconda/conda-forge only
RUN micromamba install -y -n base -c conda-forge -c bioconda bakta \
    && micromamba clean --all --yes

# Set work directory
WORKDIR /work

# Define volumes for persistent data
VOLUME ["/opt/bakta_db", "/work"]

# Default shell
CMD ["/bin/bash"]