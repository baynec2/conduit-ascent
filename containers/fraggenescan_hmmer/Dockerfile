# Base pinned by digest, not the 1.5.8 tag: a tag can be re-pushed, and this
# image's whole point is that a rebuild reproduces the same tools.
FROM mambaorg/micromamba@sha256:475730daef12ff9c0733e70092aeeefdf4c373a584c952dac3f7bdb739601990
USER root
RUN mkdir -p /work && chown -R $MAMBA_USER:$MAMBA_USER /work
USER $MAMBA_USER
# Versions pinned to what the in-use image actually contains, verified by
# inspecting /opt/conda/conda-meta in it — not chosen. An unpinned
# `micromamba install` resolves afresh on every build, so two builds of an
# identical Dockerfile produce different tools (see #68); for these packages a
# version change moves results, not just the environment.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
      fraggenescan=1.32 \
      hmmer=3.4 \
      cd-hit=4.8.1 \
      python=3.14.4 \
      biopython=1.87 \
      pandas=3.0.2 \
      plotly=6.6.0 \
      wget=1.25.0 \
      pyarrow=24.0.0 \
    && micromamba clean --all --yes
# Ensure conda base env binaries are in PATH when run via singularity exec
ENV PATH /opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
WORKDIR /work
CMD ["/bin/bash"]
