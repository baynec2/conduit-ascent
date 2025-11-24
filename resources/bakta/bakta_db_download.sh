#!/bin/bash -l
#SBATCH --job-name=bakta_db_download
#SBATCH --output=bakta_db_download_%j.out
#SBATCH --error=bakta_db_download_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=short
#SBATCH --mail-user=yac027@ucsd.edu
#SBATCH --mail-type=FAIL,END

source ~/.bashrc

source activate bakta

echo 'beginning bakta db download (full)'

bakta_db download --type full

echo 'db download finished!'

