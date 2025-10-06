#!/usr/bin/env bash
set -euo pipefail

# Simple convenience subcommand to download the MetaPhlAn database
if [[ "${1:-}" == "download-db" ]]; then
  db_dir="${2:-${METAPHLAN_DIR:-/db/metaphlan}}"
  mkdir -p "${db_dir}"
  echo "Installing MetaPhlAn database into: ${db_dir}" >&2
  micromamba run -n base metaphlan --install --bowtie2db "${db_dir}"
  exit 0
fi

# Print versions quickly
if [[ "${1:-}" == "version" || "${1:-}" == "--version" ]]; then
  exec micromamba run -n base metaphlan --version
fi

# If first argument is an option or no args provided, prepend metaphlan
if [[ ${#} -eq 0 || "${1:-}" == -* ]]; then
  set -- metaphlan "$@"
fi

exec micromamba run -n base "$@"


