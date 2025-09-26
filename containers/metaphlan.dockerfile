FROM ubuntu:22.04

LABEL base_image="ubuntu:22.04" \
      version="1" \
      software="metaphlan" \
      software.version="4.0.6" \
      description="MetaPhlAn - Metagenomic Phylogenetic Analysis for taxonomic profiling of metagenomic shotgun sequencing data" \
      homepage="https://github.com/biobakery/MetaPhlAn" \
      documentation="https://github.com/biobakery/MetaPhlAn/wiki" \
      license="https://github.com/biobakery/MetaPhlAn/blob/master/LICENSE" \
      maintainer="Charlie Bayne <baynec2@gmail.com>"

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    LANG=en_US.UTF-8 \
    LANGUAGE=en_US:en \
    LC_ALL=en_US.UTF-8 \
    PYTHONUNBUFFERED=1 \
    MPLCONFIGDIR=/tmp

# Install system dependencies and Python packages
RUN set -ex && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
        python3-dev \
        python3-venv \
        wget \
        curl \
        ca-certificates \
        locales \
        build-essential \
        libgomp1 \
        libbz2-dev \
        liblzma-dev \
        libzstd-dev \
        zlib1g-dev \
        libncurses5-dev \
        libncursesw5-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libffi-dev \
        git \
        unzip \
        gzip \
        tar && \
    # Set up locale
    locale-gen en_US.UTF-8 && \
    update-locale LANG=en_US.UTF-8 && \
    # Upgrade pip and install Python packages
    python3 -m pip install --no-cache-dir --upgrade pip setuptools wheel && \
    python3 -m pip install --no-cache-dir \
        metaphlan==4.0.6 \
        biopython \
        numpy \
        scipy \
        pandas \
        matplotlib \
        seaborn \
        plotly \
        requests \
        urllib3 && \
    # Install Bowtie2 for MetaPhlAn
    wget https://github.com/BenLangmead/bowtie2/releases/download/v2.5.2/bowtie2-2.5.2-linux-x86_64.zip && \
    unzip bowtie2-2.5.2-linux-x86_64.zip && \
    mv bowtie2-2.5.2-linux-x86_64 /opt/bowtie2 && \
    ln -s /opt/bowtie2/bowtie2 /usr/local/bin/bowtie2 && \
    ln -s /opt/bowtie2/bowtie2-build /usr/local/bin/bowtie2-build && \
    ln -s /opt/bowtie2/bowtie2-inspect /usr/local/bin/bowtie2-inspect && \
    rm bowtie2-2.5.2-linux-x86_64.zip && \
    # Install Samtools for BAM processing
    wget https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2 && \
    tar -xjf samtools-1.19.tar.bz2 && \
    cd samtools-1.19 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm -rf samtools-1.19 samtools-1.19.tar.bz2 && \
    # Clean up
    apt-get remove -y wget curl unzip gzip tar build-essential && \
    apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Create directories for MetaPhlAn databases
RUN mkdir -p /opt/metaphlan_databases && \
    chmod 755 /opt/metaphlan_databases

# Set working directory
WORKDIR /workspace

# Create a non-root user for security
RUN groupadd -r metaphlan && \
    useradd -r -g metaphlan -d /workspace -s /bin/bash metaphlan && \
    chown -R metaphlan:metaphlan /workspace /opt/metaphlan_databases

USER metaphlan

# Add MetaPhlAn to PATH
ENV PATH="/opt/bowtie2:$PATH"

# Create entrypoint script
COPY --chown=metaphlan:metaphlan <<EOF /usr/local/bin/metaphlan-entrypoint.sh
#!/bin/bash
set -e

# Default database directory
METAPHLAN_DB_DIR=\${METAPHLAN_DB_DIR:-/opt/metaphlan_databases}

# Check if database directory exists and is writable
if [ ! -d "\$METAPHLAN_DB_DIR" ]; then
    echo "Creating MetaPhlAn database directory: \$METAPHLAN_DB_DIR"
    mkdir -p "\$METAPHLAN_DB_DIR"
fi

# Set Bowtie2 database path for MetaPhlAn
export BOWTIE2_DB="\$METAPHLAN_DB_DIR"

# Execute the command
exec "\$@"
EOF

RUN chmod +x /usr/local/bin/metaphlan-entrypoint.sh

# Set entrypoint
ENTRYPOINT ["/usr/local/bin/metaphlan-entrypoint.sh"]

# Default command
CMD ["metaphlan", "--help"]

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD metaphlan --version > /dev/null || exit 1
