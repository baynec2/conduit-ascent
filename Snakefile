################################################################################
# Setup
################################################################################
# Importing necessary packages
import os
import sys
import glob
import pandas as pd
import shutil
import logging
from datetime import datetime


# Load base defaults. Any --configfile passed on the command line merges on top
# of this (later files win, per Snakemake semantics), so experiment configs only
# need to specify deltas (experiment, run_name, search_space_method, tuned flags).
configfile: "config/snakemake.yaml"


# Load configuration from command line
# Usage: snakemake --configfile path/to/config.yaml --use-apptainer
if not config.get("experiment"):
    raise ValueError("Please provide experiment name in config file, this will be used to find the correct experiment directory")

# Get experiment directory from config
EXPERIMENT_DIR = os.path.join("experiments",config["experiment"])

if not config.get("run_name"):
    raise ValueError("Please provide 'run_name' in config file to identify this analysis run")

RUN_DIR = os.path.join(EXPERIMENT_DIR, "runs", config["run_name"])

# Inject computed paths into config so all modules can read them without recomputing
config["experiment_dir"] = EXPERIMENT_DIR
config["run_dir"] = RUN_DIR

# Print the expected config file path
expected_config_path = os.path.join(RUN_DIR, "config/run_diann.cfg")

# Extracting the method from the config file
if not config.get("search_space_method"):
    raise ValueError("Please provide 'search_space_method' in config file")
# Defining the allowed methods. Will uncomment as they become supported.
ALLOWED_METHODS = [
    "ncbi_taxonomy_id",
    "uniprot_proteome_id",
    "peptidotyping",
    "unipept_hapid",
    "MAGs",
    "metaphlan",
    "hapid",
    "genome_peptidotyping",
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
print(f"Found samples: {SAMPLES}", file=sys.stderr)

# Checking to make sure that the raw file names match those in sample_annotation.txt
if not config.get("sample_annotation"):
    raise ValueError("Please provide sample_annotation file in config file")
sample_annotation = config["sample_annotation"]

# Read sample annotation file and get expected file names
try:
    sample_df = pd.read_csv(os.path.join(EXPERIMENT_DIR, sample_annotation), sep='\t')
    if 'file' not in sample_df.columns:
        raise ValueError("sample_annotation file must contain a 'file' column")
    expected_files = set(sample_df['file'].astype(str).values)
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
module setup:
  snakefile: "modules/setup/setup.smk"
  config: config
# Ways to define the search space
module metaphlan:
  snakefile: "modules/search_space/metaphlan/metaphlan.smk"
  config: config
module peptidotyping:
  snakefile: "modules/search_space/peptidotyping/peptidotyping.smk"
  config: config
module unipept_hapid:
  snakefile: "modules/search_space/unipept_hapid/unipept_hapid.smk"
  config: config
# Shared Unipept resource build (used by both peptidotyping and unipept_hapid).
module shared_unipept_resources:
  snakefile: "modules/search_space/_shared/unipept_resources.smk"
  config: config
module mags:
  snakefile: "modules/search_space/MAGs/MAGs.smk"
  config: config
module hapid:
  snakefile: "modules/search_space/hapid/hapid.smk"
  config: config
module genome_peptidotyping:
  snakefile: "modules/search_space/genome_peptidotyping/genome_peptidotyping.smk"
  config: config
# Genome download modules (optional pre-step for MAGs/hapid)
module mgnify_download:
  snakefile: "modules/genome_download/mgnify/mgnify.smk"
  config: config
module ncbi_search_space:
  snakefile: "modules/search_space/ncbi_taxonomy/ncbi_taxonomy.smk"
  config: config
module uniprot_proteome_ids_search_space:
  snakefile: "modules/search_space/uniprot_proteome_ids/uniprot_proteome_ids.smk"
  config: config
module database_processing:
  snakefile: "modules/search_space/database_processing/database_processing.smk"
  config: config
# Annotation modules
module uniprot_annotation: 
  snakefile: "modules/annotation/uniprot/annotation_uniprot.smk"
  config: config
module mag_annotation:
  snakefile: "modules/annotation/MAGs/annotation_mags.smk"
  config: config
# This provides the ability to add additional annotations from external databases.
module external_annotation:
  snakefile: "modules/annotation/external_annotations/external_annotations.smk"
  config: config
module eggnogmapper_annotation:
  snakefile: "modules/annotation/eggnogmapper/eggnogmapper.smk"
  config: config
# Diann search
module diann:
  snakefile: "modules/diann/diann.smk"
  config: config
# Building a conduit object
module build_conduit:
  snakefile: "modules/build_conduit/build_conduit.smk"
  config: config
################################################################################
# Defining all of the output files
################################################################################
rule all:
    input:
        # Reproducibility snapshot
        os.path.join(RUN_DIR, "manifest.json"),
        # Database resources
        expand(os.path.join(RUN_DIR, "database_resources/{file}"),
               file=[
                   "database.fasta",
                   #"proteome_ids.txt",
                   "taxonomy.txt",
                   "protein_info.txt",
                   "taxonomic_tree_of_database.pdf",
                   "README.md",
                   "README.html"
               ]),
        *([os.path.join(RUN_DIR, "database_resources/database.predicted.speclib")]
          if config.get("diann_search_mode", "standard") == "standard" else []),
        expand(os.path.join(RUN_DIR, "database_resources/detected_protein_resources/{file}"),
               file=[
                   "detected_protein_info.txt",
                   "detected_protein.fasta",
                   "uniprot_annotated_protein_info.txt",
                   "conduit_annotations.txt"
               ]),
        # Final output file
        conduit = os.path.join(RUN_DIR, "output_files", f"{config['experiment']}_{config['run_name']}_conduit.rds")
# Setting up the workflow. Config, apptainer, etc. 
use rule * from setup

# Search space specific workflows to generate a search space
if config["search_space_method"] == "uniprot_proteome_id":
    use rule * from uniprot_proteome_ids_search_space
    use rule * from database_processing
    use rule * from diann
    use rule * from uniprot_annotation
    use rule * from eggnogmapper_annotation
    use rule * from external_annotation

# Peptidotyping: family-level first-pass with genus fallback, then species/strain second pass
if config["search_space_method"] == "peptidotyping":
    use rule * from shared_unipept_resources
    use rule * from peptidotyping
    use rule * from ncbi_search_space
    use rule * from uniprot_proteome_ids_search_space
    use rule * from database_processing
    use rule * from diann
    use rule * from uniprot_annotation
    use rule * from eggnogmapper_annotation
    use rule * from external_annotation

# unipept_hapid: HAPiID-inspired GO-filtered first pass directly at species/strain level
if config["search_space_method"] == "unipept_hapid":
    use rule * from shared_unipept_resources
    use rule * from unipept_hapid
    use rule * from ncbi_search_space
    use rule * from uniprot_proteome_ids_search_space
    use rule * from database_processing
    use rule * from diann
    use rule * from uniprot_annotation
    use rule * from eggnogmapper_annotation
    use rule * from external_annotation

# Metaphlan feeds into the ncbi taxonomy search space
if config["search_space_method"] == "metaphlan":
    use rule * from metaphlan
    use rule * from ncbi_search_space
    use rule * from uniprot_proteome_ids_search_space
    use rule * from database_processing
    use rule * from diann
    use rule * from uniprot_annotation
    use rule * from eggnogmapper_annotation
    use rule * from external_annotation


# NCBI taxa id based workflow uses entire ncbi module.
if config["search_space_method"] == "ncbi_taxonomy_id":
    use rule * from ncbi_search_space
    use rule * from uniprot_proteome_ids_search_space
    use rule * from database_processing
    use rule * from diann
    use rule * from uniprot_annotation
    use rule * from eggnogmapper_annotation
    use rule * from external_annotation


# Genome download pre-step (runs before MAG/HAPiID if configured)
if config.get("genome_download_source") == "mgnify":
    use rule * from mgnify_download

# Search space specific workflows to generate a search space
if config["search_space_method"] == "MAGs":
    use rule * from mags
    use rule * from database_processing
    use rule * from diann
    use rule * from mag_annotation
    use rule * from eggnogmapper_annotation
    use rule * from external_annotation


# HAPiID: marker-gene profiling → greedy genome selection → MAGs DB construction
if config["search_space_method"] == "hapid":
    use rule * from hapid
    use rule * from mags
    use rule * from database_processing
    use rule * from diann
    use rule * from mag_annotation
    use rule * from eggnogmapper_annotation
    use rule * from external_annotation

# Genome peptidotyping: two-pass peptide-based detection → selected genomes → MAGs DB construction
if config["search_space_method"] == "genome_peptidotyping":
    use rule * from genome_peptidotyping
    use rule * from mags
    use rule * from database_processing
    use rule * from diann
    use rule * from mag_annotation
    use rule * from eggnogmapper_annotation
    use rule * from external_annotation


# Building Conduit Object from processed data.
use rule * from build_conduit
