################################################################################
# Setup
################################################################################
# Importing necessary packages
import os
import sys
import glob
import pandas as pd
import shlex
import shutil
import logging
from datetime import datetime
from snakemake.utils import min_version

# Pin minimum Snakemake version — this workflow uses the module system
# (introduced in 6.0) and checkpoint features that expect ≥ 8.0 semantics.
min_version("8.0")


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
    "unipept_peptidotyping",
    "unipept_hapiid",
    "genomes",
    "metaphlan",
    "hapiid",
    "genome_peptidotyping",
   # "16S"
]
# Checking that the method is allowed.
METHOD = config["search_space_method"]
# Backward-compat: the "MAGs" method was renamed to "genomes" — it accepts any
# bacterial genome FASTA (isolate assemblies, reference genomes, or MAGs), so the
# MAG-specific name was misleading. Accept the old value, warn, and normalize so
# all downstream dispatch and generated artifacts see the canonical "genomes".
if METHOD == "MAGs":
    print("WARNING: search_space_method 'MAGs' is deprecated and will be removed "
          "in a future release; use 'genomes' instead. Proceeding as 'genomes'.",
          file=sys.stderr)
    METHOD = "genomes"
    config["search_space_method"] = "genomes"
# Backward-compat: the "hapid"/"unipept_hapid" methods were misspelled; the
# canonical names are now "hapiid"/"unipept_hapiid". Accept the old spellings,
# warn, and normalize so all downstream dispatch and generated artifacts use the
# canonical name. (On-disk cache/resource paths intentionally keep the "hapid"
# spelling so existing caches are preserved.)
_HAPIID_METHOD_ALIASES = {"hapid": "hapiid", "unipept_hapid": "unipept_hapiid"}
if METHOD in _HAPIID_METHOD_ALIASES:
    _canonical = _HAPIID_METHOD_ALIASES[METHOD]
    print(f"WARNING: search_space_method '{METHOD}' is deprecated and will be "
          f"removed in a future release; use '{_canonical}' instead. "
          f"Proceeding as '{_canonical}'.", file=sys.stderr)
    METHOD = _canonical
    config["search_space_method"] = _canonical
# Backward-compat config-key aliases: the user-facing tuning keys were renamed
# hapid_* → hapiid_* (and unipept_hapid_* → unipept_hapiid_*). Internal modules
# still read the legacy key names, so populate a legacy key from its canonical
# counterpart whenever the user hasn't set the legacy name directly. This lets
# both spellings work; a user-set legacy key always wins.
_HAPIID_KEY_ALIASES = {
    "hapiid_search_mode":          "hapid_search_mode",
    "unipept_hapiid_search_mode":  "unipept_hapid_search_mode",
    "hapiid_hmm_profiles":         "hapid_hmm_profiles",
    "hapiid_percent_spectra":      "hapid_percent_spectra",
    "hapiid_infinidia_config":     "hapid_infinidia_config",
}
for _new_key, _old_key in _HAPIID_KEY_ALIASES.items():
    if _old_key not in config and _new_key in config:
        config[_old_key] = config[_new_key]
if METHOD not in ALLOWED_METHODS:
    raise ValueError(f"Method '{METHOD}' not allowed. Must be one of: {', '.join(ALLOWED_METHODS)}")

################################################################################
# DIA-NN preflight
################################################################################
# DIA-NN is supplied by the user, not shipped with this workflow — its licence
# permits one backup copy and forbids sublicensing, so it cannot live in a
# public image (issue #63). One config key, `diann_path`, names their copy;
# the helper below classifies its shape (extracted directory / bare executable
# / container image) and derives everything the rules need:
#
#   config["diann_cmd"]            -> used in shell blocks as {config[diann_cmd]}
#   config["containers"]["diann"]  -> the container (None = run on the host)
#   config["diann"]                -> provenance, lands in manifest.json
#
# Doing this at parse time means a missing or mis-shaped install fails in
# seconds with a message that says what to download and where to put it,
# rather than 40 minutes into a run. Dry runs downgrade it to a warning so the
# DAG still builds on a machine (or CI box) with no DIA-NN installed.
sys.path.insert(0, os.path.join(workflow.basedir, "modules", "_shared"))
import diann_env

# Invocations that build the DAG but never run a job. A missing DIA-NN must not
# break these: `tests/run_dry_runs.sh`, `--lint`, and the DAG-rendering step all
# run in CI, which cannot hold the binary.
_NON_EXECUTING_FLAGS = {
    "-n", "--dry-run", "--dryrun",
    "--dag", "--rulegraph", "--filegraph", "--d3dag",
    "--lint", "--list", "-l", "--list-target-rules",
    "--summary", "--detailed-summary", "--unlock",
    "--containerize", "--export-cwl", "--generate-unit-tests",
}
_IS_DRY_RUN = bool(_NON_EXECUTING_FLAGS & set(sys.argv))

try:
    DIANN = diann_env.resolve_diann_environment(config, workflow, strict=not _IS_DRY_RUN)
except diann_env.DiannEnvError as _e:
    raise ValueError(str(_e)) from None

# Publish the schema contract so module rules can pass it to the report-column
# checker without importing diann_env themselves (single writer, many readers).
config["diann_required_report_columns"] = diann_env.DIANN_REQUIRED_REPORT_COLUMNS

# strict | warn | off — see config/snakemake.yaml. Anything but "off" gates the
# cfg-snapshot rules (and therefore every DIA-NN rule) on the compatibility check.
DIANN_COMPAT_CHECK = str(config.get("diann_compat_check", "strict")).lower()
if DIANN_COMPAT_CHECK not in ("strict", "warn", "off"):
    raise ValueError(
        f"diann_compat_check must be one of strict, warn, off (got {DIANN_COMPAT_CHECK!r})"
    )
config["diann_compat_check"] = DIANN_COMPAT_CHECK
# The probe has to execute DIA-NN, so it cannot run on a dry run; skip the gate
# there rather than leaving an output that can never be produced in the DAG.
config["diann_compat_gate"] = DIANN_COMPAT_CHECK != "off" and not _IS_DRY_RUN

# Pre-render the probe invocation here rather than in the rule: building it
# needs diann_env, and the scratch dir is fixed under RUN_DIR so the command
# can be baked in (and read back in the log) instead of assembled at run time.
DIANN_PROBE_TMPDIR = os.path.join(RUN_DIR, "logs/diann/compat_probe_tmp")
config["diann_compat_probe_tmpdir"] = DIANN_PROBE_TMPDIR
config["diann_compat_probe_cmd"] = shlex.join(
    diann_env.build_compat_probe_argv(config["diann_cmd"], DIANN_PROBE_TMPDIR)
)

# Collect both .raw (Thermo native) and .mzML (open format) — DIA-NN accepts
# both as inputs. Tests ship mzML because some filtered raw files are missing
# instrument-index metadata DIA-NN needs in library-search mode; real runs
# typically use .raw directly. Directory is named ms_files/ (not raw_files/)
# because it may hold either format.
RAW_FILEPATHS = (
    glob.glob(os.path.join(EXPERIMENT_DIR, "input/ms_files/*.raw"))
    + glob.glob(os.path.join(EXPERIMENT_DIR, "input/ms_files/*.mzML"))
)
# Get just the base filenames without extension for SAMPLES
SAMPLES = []  # Initialize empty list
for filepath in RAW_FILEPATHS:
    # Get just the filename without path and without extension
    basename = os.path.basename(filepath)  # 'example2.raw' or 'example2.mzML'
    name_without_ext = os.path.splitext(basename)[0]  # 'example2'
    SAMPLES.append(name_without_ext)
# If a sample is present as both .raw and .mzML, keep a stable unique list.
SAMPLES = sorted(set(SAMPLES))

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
module unipept_peptidotyping:
  snakefile: "modules/search_space/unipept_peptidotyping/unipept_peptidotyping.smk"
  config: config
module unipept_hapid:
  snakefile: "modules/search_space/unipept_hapid/unipept_hapid.smk"
  config: config
# Shared Unipept resource build (used by both unipept_peptidotyping and unipept_hapiid).
module shared_unipept_resources:
  snakefile: "modules/search_space/_shared/unipept_resources.smk"
  config: config
module genomes:
  snakefile: "modules/search_space/genomes/genomes.smk"
  config: config
module hapid:
  snakefile: "modules/search_space/hapid/hapid.smk"
  config: config
module genome_peptidotyping:
  snakefile: "modules/search_space/genome_peptidotyping/genome_peptidotyping.smk"
  config: config
# Genome download modules (optional pre-step for genomes/hapiid)
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
module genome_annotation:
  snakefile: "modules/annotation/genomes/annotation_genomes.smk"
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
        # Database resources. NB the DIA-NN-derived intermediates (predicted
        # speclib, detected_protein_resources) are intentionally NOT listed
        # here: they are pulled transitively by build_conduit for non-empty
        # runs, and for an empty search space (no organisms detected) they must
        # NOT be required — build_conduit's checkpoint branch resolves to an
        # empty conduit without a DIA-NN search. Listing only the always-present
        # database artifacts keeps `all` satisfiable in both cases.
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

# unipept_peptidotyping: family-level first-pass with genus fallback, then species/strain second pass
if config["search_space_method"] == "unipept_peptidotyping":
    use rule * from shared_unipept_resources
    use rule * from unipept_peptidotyping
    use rule * from ncbi_search_space
    use rule * from uniprot_proteome_ids_search_space
    use rule * from database_processing
    use rule * from diann
    use rule * from uniprot_annotation
    use rule * from eggnogmapper_annotation
    use rule * from external_annotation

# unipept_hapiid: HAPiID-inspired GO-filtered first pass directly at species/strain level
if config["search_space_method"] == "unipept_hapiid":
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


# Genome download pre-step (runs before genomes/HAPiID if configured)
if config.get("genome_download_source") == "mgnify":
    use rule * from mgnify_download

# Search space specific workflows to generate a search space
if config["search_space_method"] == "genomes":
    use rule * from genomes
    use rule * from database_processing
    use rule * from diann
    use rule * from genome_annotation
    use rule * from eggnogmapper_annotation
    use rule * from external_annotation


# HAPiID: marker-gene profiling → greedy genome selection → genome DB construction
if config["search_space_method"] == "hapiid":
    use rule * from hapid
    use rule * from genomes
    use rule * from database_processing
    use rule * from diann
    use rule * from genome_annotation
    use rule * from eggnogmapper_annotation
    use rule * from external_annotation

# Genome peptidotyping: two-pass peptide-based detection → selected genomes → genome DB construction
if config["search_space_method"] == "genome_peptidotyping":
    use rule * from genome_peptidotyping
    use rule * from genomes
    use rule * from database_processing
    use rule * from diann
    use rule * from genome_annotation
    use rule * from eggnogmapper_annotation
    use rule * from external_annotation


# Building Conduit Object from processed data.
use rule * from build_conduit
