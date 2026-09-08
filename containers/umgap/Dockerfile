FROM debian:bullseye

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
    # Install Rust using rustup
    curl https://sh.rustup.rs -sSf | sh -s -- -y && \
    rustup default stable && \
    # Clone and install umgap from source
    git clone https://github.com/unipept/umgap.git /opt/umgap && \
    cd /opt/umgap && \
    cargo install --path . && \
    # Clean up
    apt-get clean && rm -rf /var/lib/apt/lists/* /opt/umgap