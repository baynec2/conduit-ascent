FROM mambaorg/micromamba:1.5.8
USER root
RUN mkdir -p /work && chown -R $MAMBA_USER:$MAMBA_USER /work
USER $MAMBA_USER
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    fraggenescan hmmer cd-hit python biopython pandas plotly wget pyarrow \
    && micromamba clean --all --yes
# Ensure conda base env binaries are in PATH when run via singularity exec
ENV PATH /opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
WORKDIR /work
CMD ["/bin/bash"]
