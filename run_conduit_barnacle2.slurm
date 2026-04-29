#!/bin/bash -l
#SBATCH --job-name=conduit
#SBATCH --error=slurm_out/conduit_%j.err
#SBATCH --output=slurm_out/conduit_%j.out
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-user=yac027@ucsd.edu
#SBATCH --mail-type=END,FAIL

# Usage:
#   sbatch run_conduit_barnacle2.slurm experiments/<exp>/config/<method>.yaml
#
# Pairs with profiles/barnacle2/ (cluster-wide snakemake flags) and
# config/snakemake.yaml (auto-loaded base defaults, per Phase 1).
set -euo pipefail

EXPERIMENT_CONFIG="${1:-}"
if [ -z "$EXPERIMENT_CONFIG" ] || [ ! -f "$EXPERIMENT_CONFIG" ]; then
    echo "ERROR: pass the experiment config YAML as the first argument." >&2
    echo "  example: sbatch $0 experiments/s77/config/hapid.yaml" >&2
    exit 1
fi

# Load and check for Singularity
module load singularity_3.6.4
singularity --version

# === Redirect Snakemake caches and temp dirs to /ddn_scratch ===
# These are user-specific (stay out of the committed profile).
export XDG_CACHE_HOME=/ddn_scratch/yac027/.cache
export TMPDIR=/ddn_scratch/yac027/tmp
export SNAKEMAKE_OUTPUT_CACHE=/ddn_scratch/yac027/.snakemake_cache

# Create directories if missing
mkdir -p $XDG_CACHE_HOME $TMPDIR $SNAKEMAKE_OUTPUT_CACHE

# Ensure barnacle2's singularity limits to match Conduit pipeline
export SINGULARITYENV_NCPUS=${SLURM_CPUS_PER_TASK:-1}
export SINGULARITYENV_OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# Run Snakemake. Profile contributes --use-singularity, --rerun-incomplete,
# and the /ddn_scratch bind. Dynamic/user-specific values are CLI args.
snakemake \
  --profile profiles/barnacle2 \
  --configfile "$EXPERIMENT_CONFIG" \
  --cores "${SLURM_CPUS_PER_TASK}" \
  --cache "${SNAKEMAKE_OUTPUT_CACHE}"

echo "FINISHED!"

