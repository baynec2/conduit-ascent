################################################################################
# Setup
################################################################################
# Importing necessary packages
import os
import glob
import pandas as pd
import shutil
import pdb  # Add this for debugging
import logging
from datetime import datetime


# Load configuration from command line
# Usage: snakemake --configfile path/to/config.yaml --sdm apptainer
if not config.get("experiment"):
    raise ValueError("Please provide experiment name in config file, this will be used to find the correct experiment directory")

# Get experiment directory from config
EXPERIMENT_DIR = os.path.join("experiments",config["experiment"])

# Print the expected config file path
expected_config_path = os.path.join(EXPERIMENT_DIR, "config/run_diann.cfg")

# Extracting the method from the config file
if not config.get("search_space_method"):
    raise ValueError("Please provide 'search_space_method' in config file")
# Defining the allowed methods. 
ALLOWED_METHODS = [
    "ncbi_taxonomy_id", # Will uncomment methods as they become supported,
   # "uniprot_proteome_id",
    "proteotyping",
   # "MAG",
    "metaphlan",
   # "16S"
]
# Checking that the method is allowed.   
METHOD = config["search_space_method"]
if METHOD not in ALLOWED_METHODS:
    raise ValueError(f"Method '{METHOD}' not allowed. Must be one of: {', '.join(ALLOWED_METHODS)}")

# Get the full paths to all raw files by joining experiment dir with pattern and then using glob
RAW_FILEPATHS = glob.glob(os.path.join(EXPERIMENT_DIR, "input/raw_files/*.raw"))
# Get just the base filenames without extension for SAMPLES
SAMPLES = []  # Initialize empty list
for filepath in RAW_FILEPATHS:
    # Get just the filename without path and without extension
    basename = os.path.basename(filepath)  # gets 'example2.raw'
    name_without_ext = os.path.splitext(basename)[0]  # gets 'example2'
    SAMPLES.append(name_without_ext)

# Print found samples for debugging
print(f"Found samples: {SAMPLES}")

# Checking to make sure that the raw file names match those in sample_annotation.txt
if not config.get("sample_annotation"):
    raise ValueError("Please provide sample_annotation file in config file")
sample_annotation = config["sample_annotation"]

# Read sample annotation file and get expected file names
try:
    sample_df = pd.read_csv(os.path.join(EXPERIMENT_DIR, sample_annotation), sep='\t')
    if 'file' not in sample_df.columns:
        raise ValueError("sample_annotation file must contain a 'file' column")
    expected_files = set(sample_df['file'].values)
except Exception as e:
    raise ValueError(f"Error reading sample annotation file: {str(e)}")

# Get actual raw file names (just the base names)
actual_files = set(SAMPLES)  # SAMPLES already contains just the base names

# Check for mismatches
missing_in_annotation = actual_files - expected_files
missing_in_raw = expected_files - actual_files

if missing_in_annotation:
    raise ValueError(f"Raw files found but not in sample annotation: {', '.join(missing_in_annotation)}")
if missing_in_raw:
    raise ValueError(f"Files in sample annotation but no matching raw files: {', '.join(missing_in_raw)}")

###############################################################################
# Module Setup and Configuration
################################################################################
if config["search_space_method"] == "ncbi_taxonomy_id":
  module search_space:
    snakefile: "modules/search_space/ncbi_taxonomy/search_space_ncbi_taxonomy.smk"
    config: config
  module annotation: 
    snakefile: "modules/annotation/ncbi_taxonomy/annotation_ncbi_taxonomy.smk"
    config: config
elif config["search_space_method"] == "proteotyping":
  module proteotyping:
    snakefile: "modules/search_space/proteotyping/generate_first_pass_proteotyping_db.smk"
    config: config
  module search_space:
    snakefile: "modules/search_space/ncbi_taxonomy/search_space_ncbi_taxonomy.smk"
    config: config
  module annotation: 
    snakefile: "modules/annotation/ncbi_taxonomy/annotation_ncbi_taxonomy.smk"
    config: config
elif config["search_space_method"] == "metaphlan":
  module search_space:
    snakefile: "modules/search_space/metaphlan/metaphlan.smk"
    config: config
  module annotation: 
    snakefile: "modules/annotation/ncbi_taxonomy/annotation_ncbi_taxonomy.smk"
    config: config
  
# elif config["search_space_method"] == "MAG":
#   module search_space:
#     snakefile: "modules/search_space/MAG/Snakefile"
#   module annotation: 
#     snakefile: "modules/annotation/MAG/Snakefile"

# elif config["search_space_method"] == "16S":
#   module search_space:
#     snakefile: "modules/search_space/16S/Snakefile"
#   module annotation: 
#     snakefile: "modules/annotation/16S/Snakefile"

# Shared modules across all methods
module setup:
  snakefile: "modules/setup/setup.smk"
  config: config
module diann:
  snakefile: "modules/diann/diann.smk"
  config: config
module build_conduit:
  snakefile: "modules/build_conduit/build_conduit.smk"
  config: config

################################################################################
# Defining all of the output files
################################################################################
rule all:
    input:
        # Config files must be created first
        os.path.join(EXPERIMENT_DIR, "config/generate_diann_spectral_library.cfg"),
        os.path.join(EXPERIMENT_DIR, "config/run_diann.cfg"),
        # Rest of the workflow outputs
        expand(os.path.join(EXPERIMENT_DIR, "output/database_resources/{file}"), 
               file=[
                   "database.fasta",
                   "proteome_ids.txt",
                   "taxonomy.txt",
                   "protein_info.txt",
                   "taxonomic_tree_of_database.pdf",
                   "database.predicted.speclib",
                   "README.md"
               ]),
        expand(os.path.join(EXPERIMENT_DIR, "output/database_resources/detected_protein_resources/{file}"),
               file=[
                   "detected_protein_info.txt",
                   "detected_protein.fasta",
                   "uniprot_annotated_protein_info.txt"
               ]),
        expand(os.path.join(EXPERIMENT_DIR, "input/raw_files/{sample}.raw"), sample=SAMPLES),       
        conduit = os.path.join(EXPERIMENT_DIR,"output","output_files",f"{config['experiment']}_conduit.rds")
# Setting up the workflow. Config, apptainer, etc. 
use rule * from setup
# Proteotyping has an additional first pass search module
if config["search_space_method"] == "proteotyping":
    use rule * from proteotyping

use rule * from search_space 
# Generating Spectral Library, Running DIA-NN, Extracting Detected Proteins
use rule * from diann
# Annotating Detected Proteins 
use rule * from annotation
# Building Conduit Object from processed data.
use rule * from build_conduit
