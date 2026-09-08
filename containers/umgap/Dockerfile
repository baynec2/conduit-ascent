# Pinned by digest, and bookworm rather than bullseye: bullseye reached
# end-of-life and its security Release file has expired, so `apt-get update`
# now fails outright ("Release file ... is expired") and this image could not
# be rebuilt at all. Digest rather than tag because `debian:bookworm` moves.
FROM debian:bookworm@sha256:6ebd97fa83deb272194a2cf015b3d26a4d538e9ad3a7a79d544c8af5b0a01443

LABEL maintainer="Charlie Bayne" \
      description="Docker container with umgap and dependencies"

# Set environment
# Pin RUSTUP_HOME/CARGO_HOME to absolute, world-readable paths so the toolchain
# is found regardless of $HOME. Under Apptainer/Singularity, $HOME is rebound to
# the calling user's host home, which would otherwise hide the default toolchain
# installed at build time and break `cargo`/`rustup` with "no default configured".
ENV RUSTUP_HOME="/usr/local/rustup" \
    CARGO_HOME="/usr/local/cargo" \
    PATH="/usr/local/cargo/bin:${PATH}"

# Pinned versions. Both are what the in-use image (umgap:e46d512, built
# 2026-06-12) actually contains, verified by inspecting it — so a rebuild
# reproduces that image rather than silently producing a different tool.
#   umgap 1.1.1 == unipept/umgap@bbb98ee4 ("Bump version number to 1.1.1")
#   rustc 1.96.0
ARG UMGAP_COMMIT=bbb98ee4
ARG RUST_VERSION=1.96.0

# Install dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl \
        git \
        uuid-runtime \
        lz4 \
        pigz \
        pv \
        build-essential \
        pkg-config \
        libssl-dev \
        ca-certificates \
        util-linux \
        unzip \
        gawk \
        libxml2-utils && \
    # Install Rust, pinned. `rustup default stable` floated the compiler: the
    # last build of this image (2026-06-12) happened to get 1.96.0, and umgap
    # upstream has since raised its minimum required Rust version, so "stable"
    # is not a fixed target.
    curl https://sh.rustup.rs -sSf | sh -s -- -y --default-toolchain ${RUST_VERSION} && \
    # Clone and install umgap at a pinned commit.
    #
    # This used to clone the default branch with no ref, so the version of the
    # LCA tool a peptidotyping run used was decided by the day the image was
    # rebuilt, and was recorded nowhere. The image in use contains umgap 1.1.1;
    # there is no v1.1.1 tag upstream (tags jump v1.0.0 -> v1.2.0), so the pin
    # is the commit that declared 1.1.1 — which is what the default branch
    # pointed at when this image was last built.
    #
    # Deliberately NOT bumped to upstream HEAD (b92e490a, umgap 1.2.0) in the
    # same commit as the base-image fix: that would ship a silent tool upgrade
    # to the peptidotyping LCA step alongside an unrelated build repair. Moving
    # to 1.2.0 is its own change, with its own revalidation.
    git clone https://github.com/unipept/umgap.git /opt/umgap && \
    cd /opt/umgap && \
    git checkout ${UMGAP_COMMIT} && \
    cargo install --path . && \
    # Clean up
    apt-get clean && rm -rf /var/lib/apt/lists/* /opt/umgap