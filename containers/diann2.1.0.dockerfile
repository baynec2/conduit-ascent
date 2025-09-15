FROM ubuntu:22.04

LABEL base_image="ubuntu:22.04" \
      version="2" \
      software="diann" \
      software.version="2.1.0" \
      description="DIA-NN - a universal software for data-independent acquisition (DIA) proteomics data processing." \
      homepage="https://github.com/vdemichev/DiaNN" \
      documentation="https://github.com/vdemichev/DiaNN" \
      license="https://github.com/vdemichev/DiaNN/LICENSE.txt" \
      maintainer="Yasset Perez-Riverol <ypriverol@gmail.com>"

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    LANG=en_US.UTF-8 \
    LANGUAGE=en_US:en \
    LC_ALL=en_US.UTF-8 \
    PATH="/usr/diann-2.1.0:$PATH" \
    DOTNET_ROOT=/usr/share/dotnet

# Install dependencies, .NET, and DIA-NN in one step
RUN set -ex && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        wget \
        unzip \
        libgomp1 \
        locales \
        ca-certificates \
        gnupg && \
    locale-gen en_US.UTF-8 && \
    update-locale LANG=en_US.UTF-8 && \
    # Install Microsoft package repo and .NET SDK 8.0
    wget https://packages.microsoft.com/config/ubuntu/22.04/packages-microsoft-prod.deb -O packages-microsoft-prod.deb && \
    dpkg -i packages-microsoft-prod.deb && \
    rm packages-microsoft-prod.deb && \
    apt-get update && \
    apt-get install -y dotnet-sdk-8.0 && \
    # Download and install DIA-NN
    wget --no-check-certificate https://github.com/vdemichev/DiaNN/releases/download/2.0/DIA-NN-2.1.0-Academia-Linux.zip -P /usr/ && \
    unzip /usr/DIA-NN-2.1.0-Academia-Linux.zip -d /usr/ && \
    rm /usr/DIA-NN-2.1.0-Academia-Linux.zip && \
    chmod +x /usr/diann-2.1.0/diann-linux && \
    ln -s /usr/diann-2.1.0/diann-linux /usr/local/bin/diann && \
    # Clean up unnecessary packages
    apt-get remove -y wget unzip && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*