# Base pinned by digest, not the 1.5.8 tag: a tag can be re-pushed, and this
# image's whole point is that a rebuild reproduces the same tools.
FROM mambaorg/micromamba@sha256:475730daef12ff9c0733e70092aeeefdf4c373a584c952dac3f7bdb739601990

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

# Install Bakta 1.12.0 from an explicit lock: every package, build and sha256
# of the environment that last ran correctly (baynec2/bakta:f5f961d). Pinning
# only the top-level tools was not enough. The rebuild that 758609a pinned kept
# bakta/diamond/blast/hmmer but re-resolved everything else: 121 of ~150
# packages changed, and perl 5.32.1 build 8 came with a relative @INC, so
# tRNAscan-SE, and with it every bakta run, died with "Can't locate strict.pm".
# Change the lock only on purpose, from an image that has been validated:
# `micromamba env export -n base --explicit`.
COPY --chown=$MAMBA_USER:$MAMBA_USER conda-linux-64.lock /tmp/conda-linux-64.lock
RUN micromamba install -y -n base --file /tmp/conda-linux-64.lock \
    && micromamba clean --all --yes

# Set work directory
WORKDIR /work

# Define volumes for persistent data
VOLUME ["/opt/bakta_db", "/work"]

# Default shell
CMD ["/bin/bash"]