#!/bin/bash -l
#SBATCH --job-name=conduit
#SBATCH --error=slurm_out/conduit_%j.err
#SBATCH --output=slurm_out/conduit_%j.out
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-user=yac027@ucsd.edu
#SBATCH --mail-type=END,FAIL

# Load and check for Singularity
module load singularity_3.6.4
singularity --version

# Ensure barnacle2's singularity limits to match Conduit pipeline
export SINGULARITYENV_NCPUS=${SLURM_CPUS_PER_TASK:-1}
export SINGULARITYENV_OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# Run Snakemake
snakemake \
  --configfile experiments/example/config/snakemake.yaml \
  --use-singularity \
  --cores ${SLURM_CPUS_PER_TASK}

echo "FINISHED!"

