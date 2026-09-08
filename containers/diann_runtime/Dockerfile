# DIA-NN *runtime* image — deliberately contains no DIA-NN.
#
# The DIA-NN licence (LICENSE.txt §3, §5) permits one backup copy and forbids
# renting/leasing/lending/sublicensing the software. Publishing an image that
# bundles the binary to a public registry is neither, so conduit-ascent no
# longer ships DIA-NN in an image (issue #63). Each user downloads DIA-NN
# themselves and accepts its terms; the workflow bind-mounts their copy into
# this image at run time.
#
# What this image provides is only the *host* environment `diann-linux` needs
# and a bare Debian/Ubuntu box does not:
#   .NET 8      the Thermo .raw readers are framework-dependent .NET assemblies
#               (see ThermoFisher.CommonCore.*.runtimeconfig.json in the zip)
#   libgomp1    OpenMP runtime that diann-linux links against
#   locales     DIA-NN parses numbers with the C/en_US locale
#
# DIA-NN's own libraries (libtorch_cpu.so, libc10.so, libtimsdata.so, ...) do
# NOT need to be here: diann-linux carries RUNPATH=$ORIGIN/, so they resolve
# from the bind-mounted directory next to the binary. No LD_LIBRARY_PATH needed.
#
# Base is ubuntu:22.04 rather than the vendor Dockerfile's debian:12 — this is
# the base the previous (DIA-NN-bundling) conduit image used and is known to
# satisfy the binary's glibc/libstdc++ requirements for DIA-NN 2.x.
FROM ubuntu:22.04

LABEL base_image="ubuntu:22.04" \
      version="1" \
      software="diann-runtime" \
      software.version="n/a" \
      description="Runtime dependencies (.NET 8, libgomp, locales) for a user-supplied DIA-NN installation. Contains no DIA-NN binary." \
      homepage="https://github.com/vdemichev/DiaNN" \
      documentation="https://github.com/baynec2/conduit-ascent" \
      maintainer="Charlie Bayne <baynec2@gmail.com>"

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=en_US.UTF-8 \
    LANGUAGE=en_US:en \
    LC_ALL=en_US.UTF-8 \
    DOTNET_ROOT=/usr/share/dotnet

RUN set -ex && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        wget \
        libgomp1 \
        locales \
        ca-certificates \
        gnupg && \
    locale-gen en_US.UTF-8 && \
    update-locale LANG=en_US.UTF-8 && \
    # Microsoft package repo + .NET SDK 8.0 (required by the Thermo .raw readers)
    wget https://packages.microsoft.com/config/ubuntu/22.04/packages-microsoft-prod.deb -O packages-microsoft-prod.deb && \
    dpkg -i packages-microsoft-prod.deb && \
    rm packages-microsoft-prod.deb && \
    apt-get update && \
    apt-get install -y dotnet-sdk-8.0 && \
    apt-get remove -y wget && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*
