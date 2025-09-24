FROM debian:bullseye

LABEL maintainer="Charlie Bayne" \
      description="Docker container with umgap and dependencies"

# Set environment
ENV PATH="/root/.cargo/bin:${PATH}"

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
    /root/.cargo/bin/rustup default stable && \
    # Clone and install umgap from source
    git clone https://github.com/unipept/umgap.git /opt/umgap && \
    cd /opt/umgap && \
    /root/.cargo/bin/cargo install --path . && \
    # Clean up
    apt-get clean && rm -rf /var/lib/apt/lists/* /opt/umgap