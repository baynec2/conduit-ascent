# Conduit: A Modular Metaproteomics Analysis Platform

Conduit is a scalable and modular workflow management system for metaproteomics data analysis, designed to seamlessly integrate with metagenomic data if available. Built using Snakemake, it provides a robust pipeline for processing Data Independent Acquisition (DIA) mass spectrometry data with particular emphasis on metaproteomics applications.

**Current version:** 0.1.0 (active development — not all features are finalized)

## Features

- **Search Space Definition**: Multiple strategies to define the protein search space for your experiment:
    - **NCBI Taxonomy IDs**: Know what taxa are in your sample? Provide NCBI taxon IDs and Conduit handles the rest.
    - **UniProt Proteome IDs**: Know specific proteome IDs? Conduit can build a search space directly from them.
    - **Peptidotyping**: A first-pass DIA-NN search using species-specific peptides to identify which taxa are present, then builds a refined search space.
    - **MetaPhlAn**: Have shotgun metagenomic data? Conduit runs MetaPhlAn profiling and uses the results to define the search space.
    - **MAGs**: Have metagenome-assembled genomes? Conduit uses Bakta to annotate them and builds the search space from the predicted proteins.
- **DIA-NN Integration**: Automated, spectral-library-free processing of DIA data.
- **Taxonomic Annotation**: Multi-level taxonomic classification of detected proteins.
- **Functional Annotation**: GO term, KEGG pathway, Pfam domain, CAZy, and eggNOG-mapper annotations.
- **R Integration**: Direct integration with R for statistical analysis and visualization via [conduitR](https://github.com/baynec2/conduitR).
- **Conduit GUI**: Explore results without writing code using [Conduit-GUI](https://github.com/baynec2/conduit-GUI).
- **Container Support**: Full containerization via Apptainer. Run on any reasonable Linux machine.
- **Named Runs**: Multiple analysis runs (e.g., different search space methods) can coexist within the same experiment directory via the `run_name` config parameter.
- **Scalable**: Runs on a single machine or scales to HPC clusters via SLURM.
- **Open Source**: MIT licensed. Customize the pipeline to fit your needs.

## Dependencies

### Required Software

- Snakemake (≥7.0.0)
- Apptainer/Singularity (≥1.1.0)

All other dependencies (R, Python, DIA-NN, MetaPhlAn, Bakta, eggNOG-mapper, etc.) are handled automatically via Apptainer containers. You only need Snakemake and Apptainer installed on your **Linux** system.

> **Note:** Conduit does not run on macOS. Windows is not recommended but may work without containers.

### Hardware Requirements

- 64 GB+ RAM recommended
- Multi-core processor (8+ cores recommended)
- Substantial storage for raw data, results, and reference databases (raw Astral MS files are ~7 GB each; reference databases for some methods can reach 50–200 GB)

## Quick Start

### 1. Install required dependencies

```bash
conda create -n conduit
conda activate conduit
conda install -c bioconda snakemake
conda install -c bioconda apptainer
```

### 2. Clone the repository

```bash
git clone https://github.com/baynec2/conduit-ascent.git
cd conduit-ascent
```

### 3. Create an experiment directory

```
experiments/your_experiment/
├── input/
│   ├── raw_files/            # Your Thermo .raw files
│   ├── ncbi_taxa_ids.txt     # NCBI taxon IDs (if using ncbi_taxonomy_id method)
│   └── sample_annotation.txt # Sample metadata
└── config/
    └── snakemake.yaml        # Experiment configuration
```

Copy `config/snakemake.yaml` as a template and fill in at minimum `experiment`, `run_name`, and `search_space_method`.

### 4. Run the workflow

```bash
snakemake --configfile experiments/your_experiment/config/snakemake.yaml --use-apptainer
```

Outputs are written to `experiments/your_experiment/runs/{run_name}/`.

### 5. Explore your results

Use `experiments/your_experiment/runs/{run_name}/output_files/{experiment}_{run_name}_conduit.rds` as input to [Conduit-GUI](https://github.com/baynec2/conduit-GUI) or [conduitR](https://github.com/baynec2/conduitR).

### 6. Repeat for additional experiments or runs

Each experiment gets its own directory. Multiple analysis runs (e.g., trying different search space methods) are namespaced via `run_name` inside the same experiment directory.

## Configuration Reference (`snakemake.yaml`)

The main configuration file controls all aspects of the workflow. Below is a full description of every parameter.

### Container Images

```yaml
containers:
  conduitr:     "docker://baynec2/conduitr:alpha"
  diann:        "docker://baynec2/diann2.1.0:alpha"
  bakta:        "docker://baynec2/bakta:alpha"
  metaphlan:    "docker://baynec2/metaphlan:alpha"
  eggnogmapper: "docker://baynec2/eggnogmapper:2.1.12"
  umgap:        "docker://baynec2/umgap:alpha"
  taxonkit:     "quay.io/biocontainers/taxonkit:0.20.0--h9ee0642_1"
```

These point to the Docker/Apptainer images used for each tool. You generally do not need to change these unless you are pinning to a specific version or using a private registry.

---

### Experiment and File Configuration

| Parameter | Description |
|-----------|-------------|
| `experiment` | **Required.** Name of the experiment. Must match the name of the directory under `experiments/`. Conduit uses this to locate all input files and write all outputs. |
| `run_name` | **Required.** Name for this specific analysis run. Outputs are written to `experiments/{experiment}/runs/{run_name}/`. Use this to run the same experiment with different methods or settings without overwriting prior results. |
| `generate_diann_spectral_library_config` | Path (relative to the main Snakefile) to the DIA-NN `.cfg` file used for spectral library generation. Defaults to `config/generate_diann_spectral_library.cfg`. |
| `run_diann_config` | Path (relative to the main Snakefile) to the DIA-NN `.cfg` file used for the main search. Defaults to `config/run_diann.cfg`. |
| `sample_annotation` | Path (relative to the experiment directory) to the sample annotation file. Defaults to `input/sample_annotation.txt`. This file maps `.raw` file names to sample metadata. |
| `output_dir` | Path (relative to the experiment directory) for the final output. Defaults to `output/`. |

---

### Search Space Configuration

| Parameter | Description |
|-----------|-------------|
| `search_space_method` | **Required.** Strategy used to define the protein search space. Options: `ncbi_taxonomy_id`, `uniprot_proteome_id`, `unipept_peptidotyping`, `MAGs`, `metaphlan`. See details below. |

#### `ncbi_taxonomy_id`
Builds the database from UniProt proteomes matching the provided NCBI taxon IDs. Requires `experiments/{experiment}/input/ncbi_taxa_ids.txt`.

#### `uniprot_proteome_id`
Builds the database from a user-supplied list of UniProt proteome IDs. Requires `experiments/{experiment}/input/proteome_ids.txt`.

#### `unipept_peptidotyping`
Performs a first-pass DIA-NN search using species-specific peptides (derived from Unipept's UMGAP LCA index) to identify which taxa are present in the sample, then builds a refined database from the detected taxa. See peptidotyping-specific parameters below.

#### `MAGs`
Accepts user-provided genome FASTAs — either metagenome-assembled genomes (MAGs) or reference genomes. Bakta annotates each FASTA and the predicted proteins form the search space. Requires `experiments/{experiment}/input/MAG_files/`.

#### `metaphlan`
Runs MetaPhlAn on shotgun metagenomic FASTQ files to profile the community, then builds a database from the detected taxa. Requires `experiments/{experiment}/input/fastq_files/`.

---

#### Peptidotyping-Specific Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `presence_min_peptides` | `3` | Minimum number of unique peptides that must be detected to call a species present. Higher values reduce false positives at the cost of sensitivity. |
| `presence_min_coverage` | `0` | Minimum protein coverage required to call a species present. Set to `0` to disable coverage filtering. |
| `min_taxon_db_peptides` | `10` | Minimum number of species-specific peptides a taxon must have in the reference database to be considered in the first-pass search. Taxa with fewer peptides are excluded as too poorly represented. |

#### MetaPhlAn-Specific Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `relative_abundance_threshold` | `0.01` | Minimum relative abundance (as a fraction, e.g., `0.01` = 1%) for a taxon to be included in the database. Increase to restrict the search space to dominant taxa; decrease to capture rarer organisms. |

---

### Universal Search Space Options

These apply to all `search_space_method` choices.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `append_additional_proteome_id` | `FALSE` | Optionally append a single extra UniProt proteome ID to the database (e.g., a host reference proteome). Common examples: Human = `UP000005640`, Mouse C57/BL6 = `UP000000589`. Set to `FALSE` to skip. |
| `append_additional_ncbi_taxa_id` | `FALSE` | Optionally append a single extra NCBI taxon ID to the database (e.g., a host species). Common examples: Human = `9606`, Mouse = `10090`. Set to `FALSE` to skip. |
| `exclude_proteome_id` | `FALSE` | Exclude a specific UniProt proteome ID from the database. Useful in peptidotyping if the first-pass incorrectly calls a contaminant. Set to `FALSE` to skip. |
| `exclude_ncbi_taxa_id` | `FALSE` | Exclude a specific NCBI taxon ID from the database. Set to `FALSE` to skip. |

---

### Resource Configuration

These parameters point to large reference databases that Conduit needs for certain methods. If the files already exist at the specified paths on your system, Conduit will use them directly. Otherwise, Conduit will download them automatically when required.

| Parameter | Default Path | Approx. Size | Description |
|-----------|-------------|-------------|-------------|
| `peptidotyping_resource_dir` | `resources/peptidotyping/` | — | Directory for peptidotyping-specific reference data. |
| `eggnog_database` | `resources/annotation/eggnog/e5.og_annotations.tsv` | ~200 MB | eggNOG ortholog annotation table used for functional annotation. |
| `pfam_db_url` | EBI FTP (current release) | — | URL for downloading the Pfam-A database. Update if a newer release is available. |
| `cazy_db_url` | dbCAN2 | — | URL for downloading the CAZy activity annotation file. |
| `eggnog_db_url` | eggNOG 5.0 | — | URL for downloading the eggNOG annotations table. |
| `metaphlan_database_dir` | `resources/metaphlan/` | ~54 GB | Directory for the MetaPhlAn reference database. Required for the `metaphlan` search space method. |
| `unipept_sequences` | `resources/databases/sequences.tsv.lz4` | ~1.46 GB | Compressed Unipept sequence table. Required for peptidotyping. |
| `unipept_taxons` | `resources/databases/taxons.tsv.lz4` | — | Compressed Unipept taxon table. Required for peptidotyping. |
| `bakta_db_dir` | `resources/bakta/db/` | ~84 GB (full) / ~2 GB (light) | Bakta annotation database. Required for the `MAGs` method. Controlled by `bakta_db_type`. |
| `bakta_db_type` | `"light"` | — | Which Bakta database to use: `"full"` (~84 GB) or `"light"` (~2 GB). The light database is faster to download but less comprehensive. |
| `eggnogmapper_db_dir` | `resources/eggnogmapper/db/` | ~50 GB (full) / ~8 GB (bacteria) | eggNOG-mapper database directory for sequence-based functional annotation. |
| `gbtk_database` | `resources/gbtk/` | ~140 GB | GTDB-Tk database. Reserved for future MAG taxonomy use. |

---

## Project Structure

```
conduit-ascent/
├── Snakefile                         # Main workflow orchestration
├── VERSION                           # Current version
├── config/                           # Default configuration files (templates)
│   ├── snakemake.yaml                # Main config template
│   ├── generate_diann_spectral_library.cfg
│   ├── run_diann.cfg
│   ├── peptidotyping_firstpass_diann.cfg
│   └── peptidotyping_firstpass_diann_spectral_library.cfg
├── modules/                          # Snakemake modules
│   ├── setup/                        # Workflow setup (config copying, image building)
│   ├── search_space/                 # Search space definition
│   │   ├── ncbi_taxonomy/            # From NCBI taxonomy IDs
│   │   ├── uniprot_proteome_ids/     # From UniProt proteome IDs
│   │   ├── peptidotyping/            # First-pass species detection
│   │   ├── metaphlan/                # From MetaPhlAn metagenomic profiling
│   │   ├── MAGs/                     # From metagenome-assembled genomes
│   │   └── database_processing/      # Shared post-processing
│   ├── diann/                        # DIA-NN identification and quantification
│   ├── annotation/                   # Protein and taxonomic annotation
│   │   ├── uniprot/                  # UniProt-based annotation
│   │   ├── MAGs/                     # Bakta + UniProt annotation for MAG proteins
│   │   ├── eggnogmapper/             # eggNOG-mapper functional annotation
│   │   └── external_annotations/     # KEGG, Pfam, CAZy annotations
│   └── build_conduit/                # Builds final Conduit RDS object
├── containers/                       # Apptainer/Docker container definitions
│   ├── bakta/
│   ├── conduitR/
│   ├── diann/
│   ├── eggnogmapper/
│   ├── metaphlan/
│   └── umgap/
├── experiments/                      # Experiment directories
│   └── {experiment_name}/
│       ├── config/                   # Per-experiment config
│       ├── input/
│       │   ├── raw_files/            # Thermo .raw MS files
│       │   ├── sample_annotation.txt # Sample metadata
│       │   ├── ncbi_taxa_ids.txt     # Taxon IDs (ncbi_taxonomy_id method)
│       │   ├── proteome_ids.txt      # Proteome IDs (uniprot_proteome_id method)
│       │   ├── fastq_files/          # FASTQ files (metaphlan method)
│       │   └── MAG_files/            # MAG FASTA files (MAGs method)
│       └── runs/
│           └── {run_name}/           # All outputs for a given run
│               ├── config/
│               ├── database_resources/
│               ├── logs/
│               └── output_files/
├── tests/                            # Test suite
│   └── configs/                      # Integration test configs
└── .github/workflows/                # CI (dry-run tests)
```

## Output Structure

All outputs for a run are written to `experiments/{experiment}/runs/{run_name}/`:

| Output | Description |
|--------|-------------|
| `database_resources/database.fasta` | Protein FASTA used as the search space |
| `database_resources/taxonomy.txt` | Taxonomy information for all organisms in the database |
| `database_resources/protein_info.txt` | Per-protein metadata |
| `database_resources/proteome_ids.txt` | UniProt proteome IDs included in the database |
| `database_resources/database.predicted.speclib` | DIA-NN predicted spectral library |
| `database_resources/detected_protein_resources/` | Proteins detected by DIA-NN and their annotations |
| `output_files/{experiment}_{run_name}_conduit.rds` | Final Conduit R object for use with conduitR or Conduit-GUI |

## Running Conduit on Barnacle2 (Knight Lab HPC)

These instructions are for Knight Lab members running Conduit on Barnacle2 via SLURM. They can be adapted to other SLURM-based HPC systems.

> For tutorial purposes you will also need files in `experiments/example/input/database_resources` and `experiments/example/input/raw_files`. Contact baynec2 directly for these files.

### 1. Login

```bash
ssh <username>@barnacle2.ucsd.edu
```

### 2. Clone the repository onto scratch

Clone into your `/ddn_scratch` space — it's large and not purged, so the repo,
your MS files, all outputs (`runs/`), and the (re)built peptidotyping resources
live there and stay visible inside the rule containers (the profile binds
`/ddn_scratch`). Avoid `$HOME`, which is quota-limited.

```bash
cd /ddn_scratch/$USER
git clone https://github.com/baynec2/conduit-ascent.git
cd conduit-ascent
```

### 3. Install Miniforge and create the Snakemake environment

Barnacle2 has **no conda/Miniforge module** (`module avail` lists only tools like
`singularity_3.6.4`), so install Miniforge yourself, into `$HOME` — a Miniforge
install plus this env is small (a few hundred MB), and `$HOME` is backed up and
mounted on the compute nodes. Miniforge defaults to the conda-forge channel and
ships `mamba` as a fast solver (`mamba` and `conda` are interchangeable below).

```bash
cd ~
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b            # installs to $HOME/miniforge3
source $HOME/miniforge3/etc/profile.d/conda.sh   # activate in this shell
$HOME/miniforge3/bin/conda init bash             # load on future logins
mamba --version                                  # sanity check
```

Then create a dedicated environment with Snakemake and the SLURM executor plugin.
Don't install singularity/apptainer into it — Barnacle2 provides Singularity as a
module (loaded per session in step 4), and keeping it out of the env leaves the
module's Singularity on `PATH`.

```bash
mamba create -n conduit -c conda-forge -c bioconda \
    snakemake snakemake-executor-plugin-slurm
conda activate conduit
snakemake --version                      # sanity check
```

The `snakemake-executor-plugin-slurm` package is what lets Snakemake submit each
pipeline rule as its own SLURM job (see `profiles/barnacle2/`), rather than
running everything inside one large allocation.

### 4. Run the workflow

The `profiles/barnacle2/` profile uses the SLURM **executor**: Snakemake itself
submits each rule as its own SLURM job. The Snakemake process is lightweight (it
just submits and polls jobs), so run it directly from a login node inside a
`tmux`/`screen` session so it survives your SSH disconnecting.

```bash
tmux new -s conduit          # so the run survives disconnects

# --- one-time-per-session setup ---
conda activate conduit
module load singularity_3.6.4
singularity --version       # confirm Singularity resolves to the module

# Point Snakemake's caches/temp at your scratch space. SLURM exports this
# environment (incl. the loaded module) to the per-rule jobs, so Singularity
# is on PATH inside them too.
export XDG_CACHE_HOME="/ddn_scratch/${USER}/.cache"
export TMPDIR="/ddn_scratch/${USER}/tmp"
export SNAKEMAKE_OUTPUT_CACHE="/ddn_scratch/${USER}/.snakemake_cache"
mkdir -p "$XDG_CACHE_HOME" "$TMPDIR" "$SNAKEMAKE_OUTPUT_CACHE"

# --- launch the workflow ---
snakemake \
  --profile profiles/barnacle2 \
  --configfile experiments/<exp>/config/<method>.yaml \
  --cache "$SNAKEMAKE_OUTPUT_CACHE"
```

Detach from `tmux` with `Ctrl-b d`; reattach later with `tmux attach -t conduit`.
Check the per-rule SLURM jobs Snakemake has submitted with `squeue --me`, and
per-rule logs under `runs/{run_name}/logs/` (workflow) and
`.snakemake/slurm_logs/` (raw SLURM stdout/stderr, written by the executor).

The profile caps concurrent SLURM jobs at 50, lets each job use up to a full
node (64 cores), and scales memory at ~8 GB/core (capped ~24 GB below the node's
514 GB). Walltime is left to barnacle2's partition default — `short` allows 14
days, which covers even the multi-day peptidotyping resource rebuild. After a
run, the `benchmarks/*.tsv` files record each rule's real peak memory and
runtime; use them to tighten the `set-resources` values in the profile.

## Troubleshooting

### NCBI API failures
If the `ncbi_taxonomy_id` method fails to fetch taxonomy, the NCBI API may be temporarily unavailable. Wait a few minutes and retry.

### Configuration errors
Ensure your `snakemake.yaml` has values for all required fields: `experiment`, `run_name`, `search_space_method`, and `sample_annotation`.

### Getting help
- Check the [Issues](https://github.com/baynec2/conduit-ascent/issues) page for known problems
- Open a new issue with your error message, configuration file, and relevant log from `runs/{run_name}/logs/`

## Misc

### Visualize the DAG

```bash
snakemake --configfile experiments/example/config/snakemake.yaml --dag \
  | grep -v '^Found samples' \
  | dot -Tpdf > dag.pdf
```

## Contributing

Contributions are welcome! Please open an issue or pull request.

### New Module Development Guide

#### Module Structure

- Each module lives in its own directory describing its task (e.g., `search_space/`, `diann/`, `annotation/`).
- Alternative implementations of the same task go in subdirectories under the task directory.
- Each module must contain a Snakemake file named `{module_name}.smk` and a `scripts/` subdirectory.
- Logs must be written to the main experiment's `logs/` directory.
- `rule all` must only be defined in the top-level `Snakefile`, never in modules.

#### Search Space and Annotation Pairing

Each search space module **must have a corresponding annotation module**, as different search space strategies require different protein annotation approaches.

#### Required Outputs: Search Space Modules

All outputs must be placed in `database_resources/`. If a file cannot be generated for a given method, produce it with `NA` values.

| File | Contents |
|------|----------|
| `database.fasta` | All protein sequences in the search space, UniProt-style headers |
| `proteome_ids.txt` | UniProt proteome IDs used |
| `taxonomy.txt` | Taxonomy information for the database |
| `protein_info.txt` | Per-protein metadata |
| `taxonomic_tree_of_database.pdf` | Taxonomic tree visualization of the search space |
| `database.predicted.speclib` | DIA-NN predicted spectral library |
| `README.md` | Metrics about the database resources |

#### Required Outputs: Annotation Modules

| File | Contents |
|------|----------|
| `detected_protein_info.txt` | `protein_info.txt` filtered to DIA-NN-detected proteins |
| `detected_protein.fasta` | Detected proteins in FASTA format |
| `annotated_protein_info.txt` | Detected proteins with annotations added |
| `go_annotations.txt` | GO term annotations in long format |
| `subcellular_locations.txt` | Subcellular location predictions |
| `kegg_annotations.txt` | KEGG pathway annotations |

## License

This project is licensed under the MIT License — see the LICENSE file for details.

## Acknowledgments

Conduit would not be possible without the great work of many people.

- Snakemake developers
- DIA-NN developers
- R Bioconductor community
- R tidyverse community
- UniProt consortium
- NCBI taxonomy database maintainers
- KEGG database maintainers
- Gene Ontology consortium
- MetaPhlAn team
- Bakta developers
- eggNOG-mapper team
- Unipept team
