# Conduit: A Modular Metaproteomics Analysis Platform

Conduit is a scalable and modular workflow management system for metaproteomics data analysis, designed to seamlessly integrate with metagenomic data if available. Built using Snakemake, it provides a robust pipeline for processing Data Independent Acquisition (DIA) mass spectrometry data with particular emphasis on metaproteomics applications.

**Current version:** 0.1.0 (active development — not all features are finalized)

## Features

- **Search Space Definition**: Multiple strategies to define the protein search space for your experiment:
    - **NCBI Taxonomy IDs**: Know what taxa are in your sample? Provide NCBI taxon IDs and Conduit handles the rest.
    - **UniProt Proteome IDs**: Know specific proteome IDs? Conduit can build a search space directly from them.
    - **Peptidotyping**: A first-pass DIA-NN search using species-specific peptides to identify which taxa are present, then builds a refined search space. Available against UniProt-derived peptides (`unipept_peptidotyping`), a HAPiID-style species-first variant (`unipept_hapiid`), or peptides digested from your own genomes with no UniProt lookup (`genome_peptidotyping`).
    - **MetaPhlAn**: Have shotgun metagenomic data? Conduit runs MetaPhlAn profiling and uses the results to define the search space.
    - **Genomes**: Have bacterial genome FASTAs (metagenome-assembled or reference)? Conduit uses Bakta to annotate them and builds the search space from the predicted proteins. Genomes can also be downloaded from an MGnify catalog instead of supplied by hand.
    - **HAPiID**: Marker-gene profiling across your genomes, then greedy selection of the smallest genome set covering most of the annotated spectra.
- **DIA-NN Integration**: Automated, spectral-library-free processing of DIA data.
- **Taxonomic Annotation**: Multi-level taxonomic classification of detected proteins.
- **Functional Annotation**: GO term, KEGG pathway, Pfam domain, CAZy, and eggNOG-mapper annotations.
- **R Integration**: Direct integration with R for statistical analysis and visualization via [conduitR](https://github.com/baynec2/conduitR).
- **Conduit GUI**: Explore results without writing code using [Conduit-GUI](https://github.com/baynec2/conduit-GUI).
- **Container Support**: Full containerization via Apptainer. Run on any reasonable Linux machine.
- **Named Runs**: Multiple analysis runs (e.g., different search space methods) can coexist within the same experiment directory via the `run_name` config parameter.
- **Scalable**: Runs on a single machine or scales to HPC clusters via SLURM.
- **Readable and modular**: every stage is a Snakemake rule you can read and adapt. Released under the [MIT License](LICENSE).

## Dependencies

### Required Software

- Snakemake (≥7.0.0)
- Apptainer/Singularity (≥1.1.0)
- **DIA-NN** — you download this yourself; see [Obtaining DIA-NN](#obtaining-dia-nn)

All other dependencies (R, Python, MetaPhlAn, Bakta, eggNOG-mapper, etc.) are handled automatically via Apptainer containers. You only need Snakemake, Apptainer, and your own copy of DIA-NN on your **Linux** system.

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

### 2b. Obtain DIA-NN

Conduit does not ship DIA-NN. Download it, unzip it into `resources/diann/`, and you
are done — see [Obtaining DIA-NN](#obtaining-dia-nn) for the details, including why the
download URL says `2.0` no matter which version you are fetching.

```bash
mkdir -p resources/diann
curl -fL -o /tmp/diann.zip \
  https://github.com/vdemichev/DiaNN/releases/download/2.0/DIA-NN-2.5.0-Academia-Linux.zip
unzip /tmp/diann.zip -d resources/diann && rm /tmp/diann.zip
chmod +x resources/diann/diann-2.5.0/diann-linux
```

By downloading DIA-NN you are accepting its licence directly from its authors; free for
academic use, commercial use requires a licence from them.

### 3. Create an experiment directory

```
experiments/your_experiment/
├── input/
│   ├── ms_files/             # Your .raw (Thermo) or .mzML files
│   ├── ncbi_taxa_ids.txt     # NCBI taxon IDs (if using ncbi_taxonomy_id method)
│   └── sample_annotation.txt # Sample metadata; needs a `file` column naming
│                             # each MS file without its extension
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
  conduitr:           "docker://baynec2/conduitr:a5cfeae"
  diann_runtime:      "docker://baynec2/diann_runtime:<sha>"
  bakta:              "docker://baynec2/bakta:f5f961d"
  metaphlan:          "docker://baynec2/metaphlan:f5f961d"
  eggnogmapper:       "docker://baynec2/eggnogmapper:f5f961d"
  umgap:              "docker://baynec2/umgap:e46d512"
  fraggenescan_hmmer: "docker://baynec2/fraggenescan_hmmer:f5f961d"
  taxonkit:           "docker://quay.io/biocontainers/taxonkit:0.20.0--h9ee0642_1"
```

These point to the Docker/Apptainer images used for each tool. You generally do not need to
change them unless you are pinning to a specific version or using a private registry.

There is deliberately no `diann` key. DIA-NN is supplied by you (see
[Obtaining DIA-NN](#obtaining-dia-nn)) and the workflow derives its container from
`diann_path` at parse time. `diann_runtime` holds only DIA-NN's *host* dependencies
(.NET 8, libgomp, locales) — no DIA-NN binary — and is what your extracted DIA-NN
directory gets bind-mounted into.

Tags are short commit SHAs rather than moving tags like `:latest`. Apptainer caches images
by URI, so a moving tag would not refresh once cached and you would silently keep running
an old image; a SHA makes the change explicit. See the tag-bumping note in `CLAUDE.md`
before starting a large run during active development.

---

### Experiment and File Configuration

| Parameter | Description |
|-----------|-------------|
| `experiment` | **Required.** Name of the experiment. Must match the name of the directory under `experiments/`. Conduit uses this to locate all input files and write all outputs. |
| `run_name` | **Required.** Name for this specific analysis run. Outputs are written to `experiments/{experiment}/runs/{run_name}/`. Use this to run the same experiment with different methods or settings without overwriting prior results. |
| `sample_annotation` | Path (relative to the experiment directory) to the sample annotation file. Defaults to `input/sample_annotation.txt`. Maps MS file names (without extension) to sample metadata via a required `file` column. **Effectively fixed:** `build_conduit.smk` hardcodes the default path and ignores this key, so overriding it passes the Snakefile's checks and then fails downstream. |
| `diann_search_mode` | `standard` (default) or `infinidia`. `standard` runs a three-stage library / per-file / combine search; `infinidia` runs a single monolithic search with `--pre-search --pre-filter`. Also decides whether a predicted spectral library is built at all. |

#### DIA-NN configuration files

Each key points at a `.cfg` file of DIA-NN flags. All of them are snapshotted into
`runs/{run_name}/config/` at run start, so the settings a run used are frozen with its
outputs. Override a key in your experiment YAML to swap in a variant without editing the
shared file.

| Parameter | Default | Used by |
|-----------|---------|---------|
| `run_diann_config` | `config/run_diann.cfg` | the main search, every method (when `diann_search_mode: standard`) |
| `diann_spectral_library_base_config` | `config/diann_spectral_library_base.cfg` | every spectral-library prediction step |
| `diann_library_search_base_config` | `config/diann_library_search_base.cfg` | `hapiid`, `unipept_hapiid` in standard mode |
| `hapiid_infinidia_config` | `config/hapid_infinidia.cfg` | `hapiid`, `unipept_hapiid` in InfiniDIA mode |
| `peptidotyping_standard_config` | `config/peptidotyping_standard.cfg` | `unipept_peptidotyping`, `genome_peptidotyping` in standard mode |
| `peptidotyping_infinidia_config` | `config/peptidotyping_infinidia.cfg` | `unipept_peptidotyping`, `genome_peptidotyping` in InfiniDIA mode |

Note that inputs, outputs, scratch paths and thread counts (`--fasta`, `--dir`, `--lib`,
`--out`, `--out-lib`, `--temp`, `--threads`) are supplied by the rules themselves, and
several rules also pass digest flags inline. Setting those in a `.cfg` has no effect.

---

### Search Space Configuration

| Parameter | Description |
|-----------|-------------|
| `search_space_method` | **Required.** Strategy used to define the protein search space. Options: `ncbi_taxonomy_id`, `uniprot_proteome_id`, `unipept_peptidotyping`, `unipept_hapiid`, `genomes`, `metaphlan`, `hapiid`, `genome_peptidotyping`. See details below. |

Every method requires `experiments/{experiment}/input/ms_files/` (`.raw` or `.mzML`) and
`experiments/{experiment}/input/sample_annotation.txt`. The requirements listed per method
below are in addition to those.

#### `ncbi_taxonomy_id`
Builds the database from UniProt proteomes matching the provided NCBI taxon IDs. Requires `experiments/{experiment}/input/ncbi_taxa_ids.txt`.

#### `uniprot_proteome_id`
Builds the database from a user-supplied list of UniProt proteome IDs. Requires `experiments/{experiment}/input/proteome_ids.txt`.

#### `unipept_peptidotyping`
Performs a first-pass DIA-NN search using species-specific peptides (derived from Unipept's UMGAP LCA index) to identify which taxa are present in the sample, then builds a refined database from the detected taxa. See peptidotyping-specific parameters below.

#### `unipept_hapiid`
HAPiID-inspired variant of the above: a GO-filtered first pass straight at species/strain level, rather than family-first with a species second pass. The old name `unipept_hapid` is still accepted as a deprecated alias.

#### `genomes`
Accepts user-provided genome FASTAs — either metagenome-assembled genomes (MAGs) or reference genomes. Bakta annotates each FASTA and the predicted proteins form the search space. The old method name `MAGs` is still accepted as a deprecated alias.

Requires `experiments/{experiment}/input/genome_files/` containing:

- one or more `.fa`, `.fna` or `.fasta` genome files, and
- **`taxonomy.txt`** — a tab-separated table with a required `genome` column matching the
  FASTA base names, plus optional `domain kingdom phylum class order family genus species`
  columns (missing ranks default to `NA`).

Omitting `taxonomy.txt` is the most common first-run failure. The legacy `MAG_files/`
directory name is still accepted when `genome_files/` is absent.

Both requirements are lifted when `genome_download_source: mgnify` — the catalog supplies
the genomes and the taxonomy table instead.

#### `hapiid`
Marker-gene (ribosomal protein / elongation factor) profiling across your genomes, then greedy selection of the smallest genome set covering `hapiid_percent_spectra` of the annotated spectra. Takes the same `genome_files/` inputs as `genomes`, and additionally needs the HMM profiles at `hapiid_hmm_profiles`. The old name `hapid` is still accepted as a deprecated alias.

#### `genome_peptidotyping`
Two-pass peptide-based detection over your genomes with no UniProt lookup: tryptic peptides are digested from the genomes, assigned an LCA, and searched to decide which genomes are present. Takes the same `genome_files/` inputs as `genomes`. The rank columns in `taxonomy.txt` are functionally required here — the LCA computation and the species-to-genome match both read them.

#### `metaphlan`
Runs MetaPhlAn on shotgun metagenomic FASTQ files to profile the community, then builds a database from the detected taxa. Requires `experiments/{experiment}/input/fastq_files/`.

**Reads must be gzipped.** The rule globs `*.fastq.gz` only; an uncompressed `.fastq` is
not matched, so it contributes no sample and the run continues with an empty profile and
no error.

---

#### Peptidotyping-Specific Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
Applies to `unipept_peptidotyping` and `genome_peptidotyping`. Both methods run two
detection passes, and every key below exists once per pass — the first-pass name is given
here, and the second-pass twin is the same name with `second` in place of `first` (e.g.
`peptidotyping_second_pass_method`).

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_taxon_db_peptides` | `10` | Minimum number of species-specific peptides a taxon must have in the reference database to enter the detection search at all. Taxa with fewer are excluded as too poorly represented. Not per-pass. |
| `peptidotyping_first_pass_method` | `count` | The presence-call rule. `count` is recommended: a taxon is present if it has at least `..._min_confident_peptides` distinct peptides at DIA-NN `Q.Value <= 0.01`, resting on DIA-NN's own validated precursor FDR. `enrichment` and `qvalue` both depend on the reported decoy null, which InfiniDIA's `--pre-filter` censors, and are kept for provenance rather than use. |
| `peptidotyping_first_pass_min_confident_peptides` | `10` | Under `method: count`, how many confident peptides a taxon needs. Raising it tightens the "congener halo" of close relatives at some cost to low-abundance members. |
| `peptidotyping_first_pass_min_peptides` | `10` | Minimum distinct contributing peptides regardless of confidence. Raises the bar without regard to quality. |
| `peptidotyping_first_pass_margin` | `2.0` | Under `method: enrichment` only, the multiple of the decoy noise rate a lineage's per-peptide score must clear. Lowering it enlarges the second-pass database, which is searched library-free and can exhaust memory. |
| `peptidotyping_first_pass_qvalue_threshold` | `0.05` | Under `method: qvalue` only, the picked q-value cutoff. |
| `peptidotyping_first_pass_score_fraction_threshold` | `null` | Optional score-coverage gate: keep the smallest top-by-score prefix of taxa reaching this fraction of total score. Disabled by default because it drops real low-abundance taxa. |
| `peptidotyping_first_pass_max_taxa` | `null` | Optional cap: keep the top N taxa by score. **Set at most one** of this and `..._score_fraction_threshold` for the same pass — setting both is an error. |

#### Per-method search modes

Independent of `diann_search_mode`, which controls only the main search. Each is
`standard` (build a predicted spectral library, then search with `--lib`) or `infinidia`
(search the FASTA directly with `--pre-search --pre-filter`).

| Parameter | Default |
|-----------|---------|
| `unipept_peptidotyping_search_mode` | `infinidia` |
| `genome_peptidotyping_search_mode` | `infinidia` |
| `unipept_hapiid_search_mode` | `standard` |
| `hapiid_search_mode` | `standard` |

#### HAPiID-Specific Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `hapiid_percent_spectra` | `80` | Take the smallest set of genomes covering this percentage of marker-gene-annotated spectra. Applies to `hapiid` and `unipept_hapiid`. |
| `hapiid_hmm_profiles` | `resources/hapid/ribP_elonF_profiles_refined_manually.hmm` | HMM profiles for marker-gene identification. |

#### Genome Source (optional MGnify download)

Applies to `genomes`, `hapiid` and `genome_peptidotyping`.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `genome_download_source` | `FALSE` | `FALSE` uses the genomes you provide. `mgnify` downloads a catalog instead, replacing **both** your FASTAs and your `taxonomy.txt`. |
| `mgnify_catalog` | `""` | Catalog and version, e.g. `human-gut/v2.0.2`. The string is the version pin. Required when the source is `mgnify`; it is not validated, and a typo surfaces late as a parse failure on the downloaded metadata. |
| `mgnify_taxonomy_filter` | `FALSE` | Optional GTDB-lineage substring filter, e.g. `p__Firmicutes`. |
| `mgnify_max_genomes` | `0` | Cap on species representatives to download. `0` means no limit. |
| `mgnify_ftp_base` | EBI FTP | Base URL; should not need changing. |

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
| `append_additional_ncbi_taxa_id` | `FALSE` | Optionally append a single extra NCBI taxon ID to the database (e.g., a host species). Common examples: Human = `9606`, Mouse = `10090`. Set to `FALSE` to skip. Has no effect under `uniprot_proteome_id`, which never routes through the NCBI taxonomy module. |

On the genome-sourced methods (`genomes`, `hapiid`, `genome_peptidotyping`) the two keys
are **mutually exclusive** — setting both raises. On the UniProt-based methods they apply
at different stages and stack.

---

### Resource Configuration

These parameters point to large reference databases that Conduit needs for certain methods. If the files already exist at the specified paths on your system, Conduit will use them directly. Otherwise, Conduit will download them automatically when required.

| Parameter | Default Path | Approx. Size | Description |
|-----------|-------------|-------------|-------------|
| `peptidotyping_resource_dir` | `resources/peptidotyping/` | ~170 GB | UMGAP-derived `sequences.tsv.lz4` + `taxons.tsv.lz4`. Required for `unipept_peptidotyping` and `unipept_hapiid`. |
| `taxonkit_db_dir` | `resources/peptidotyping/taxonkit` | — | NCBI taxonomy database for taxonkit. Conventionally lives under `peptidotyping_resource_dir`. |
| `metaphlan_database_dir` | `resources/metaphlan/` | ~54 GB | MetaPhlAn reference database. Required for the `metaphlan` method. |
| `bakta_db_dir` | `resources/bakta/db/` | ~84 GB (full) / ~2 GB (light) | Bakta annotation database. Required for the genome-sourced methods. Controlled by `bakta_db_type`. |
| `bakta_db_type` | `"light"` | — | Which Bakta database to use: `"full"` (~84 GB) or `"light"` (~2 GB). Light is sufficient for most bacterial work; use full only when you need the plasmid/viral databases. |
| `eggnogmapper_db_dir` | `resources/eggnogmapper/db/` | ~50 GB (full) / ~8 GB (bacteria) | eggNOG-mapper database directory for sequence-based functional annotation. |
| `genome_resource_dir` | `resources/genome_databases/derived` | — | Cache of per-genome artifacts (Prodigal/FGS, HMMER, bakta). Shared across experiments only when the source is MGnify, whose accessions are globally unique. |
| `genome_set_resource_dir` | `resources/genome_databases/derived_sets` | — | Cache of per-genome-set aggregates (LCA peptide databases, HAPiID libraries). |
| `mgnify_cache_dir` | `resources/genome_databases/mgnify` | — | Downloaded MGnify catalogs, one subdirectory per catalog. |
| `pfam_db_url` | EBI FTP (current release) | — | URL for downloading the Pfam-A database. |
| `cazy_db_url` | dbCAN2 | — | URL for downloading the CAZy activity annotation file. |
| `eggnog_db_url` | eggNOG 5.0 | — | URL for downloading the eggNOG annotations table. |
| `go_obo_url`, `kegg_rest_url`, `enzyme_dat_url`, `pfam_clans_url` | pinned releases | — | Term-name dictionaries used to fill annotation descriptions. |

These paths are machine-specific, so they belong in your Snakemake **profile** rather than
in an experiment config — the profile's `config:` block overrides both the base config and
the experiment config. See `profiles/nanopore-catalyst/config.yaml` for a worked example.

---

## Obtaining DIA-NN

Conduit does not distribute DIA-NN, and no Conduit container contains it.

DIA-NN's licence permits **one copy for backup purposes** and forbids renting, leasing,
lending, or sublicensing the software. Publishing a public registry image with the binary
inside is none of those things, so the binary is not ours to ship. You download it yourself
and accept its terms directly — which is where licence acceptance belongs anyway. DIA-NN is
free for academic use; commercial use requires a licence from its authors.

- Download: <https://github.com/vdemichev/DiaNN/releases> (`DIA-NN-<version>-Academia-Linux.zip`)
- Licence: <https://github.com/vdemichev/DiaNN/blob/master/LICENSE.txt>

> **The releases page is misleading — read this before you download.** The newest *release*
> shown is **2.0**, dated January 2025. That is not the newest DIA-NN. Every build since —
> 2.0.1 through 2.6.1, twelve Linux builds at the time of writing — is published as an
> **asset attached to that same `2.0` tag**, not as its own release. So the download URL
> always carries `2.0` in the path regardless of which version you are fetching:
>
> ```
> https://github.com/vdemichev/DiaNN/releases/download/2.0/DIA-NN-2.5.0-Academia-Linux.zip
>                                             ^^^ always 2.0, never the version you want
> ```
>
> Scroll to the **Assets** list on the `2.0` release to see what is actually available, or
> list them from the command line:
>
> ```bash
> gh api repos/vdemichev/DiaNN/releases/tags/2.0 \
>   --jq '.assets[].name | select(test("Academia-Linux"))'
> ```
>
> Taking the release title at face value lands you on 2.0, which Conduit rejects as below
> the 2.2 minimum.

Tested versions: **2.3.0** and **2.5.0**. Anything below **2.2** is rejected — earlier
builds silently ignore `--pre-search` / `--pre-filter`, which would turn every InfiniDIA
search into a plain library search with no error. See
[What gets checked, and when](#what-gets-checked-and-when) for what "tested" does and does
not mean.

### Where to put it

One config key, `diann_path`, names your copy. Three shapes are autodetected:

| Shape | Detected when | What happens |
|---|---|---|
| **A — directory** *(recommended)* | the path is a directory containing `diann-linux` | Runs inside the `diann_runtime` image with your directory bind-mounted in |
| **B — executable** | the path is an executable file | Runs directly on the host, no container. Use for `DIA-NN.AppImage` (self-contained) or a host that already has .NET 8 |
| **C — image** | the path starts with `docker://` or ends in `.sif` | Used as the container as-is. Migration path if you already built your own full image |

Mode A is what most people want:

```bash
mkdir -p resources/diann
unzip DIA-NN-2.5.0-Academia-Linux.zip -d resources/diann
chmod +x resources/diann/diann-2.5.0/diann-linux
```

which gives the default layout, so nothing else needs configuring:

```
resources/diann/diann-2.5.0/
├── diann-linux                      # the CLI Conduit invokes
├── libtorch_cpu.so, libc10.so, ...  # bundled, found via RUNPATH=$ORIGIN
└── models/
```

`diann-linux` is a framework-dependent .NET binary and also needs libgomp, so pointing at
it on a bare host is not enough on its own. That is what `containers.diann_runtime`
supplies: Debian/Ubuntu + .NET 8 + libgomp + locales, and **no DIA-NN**. Your directory is
bind-mounted in and the binary is executed from the mount, so nothing license-restricted
ever enters an image we publish — and you avoid an `apptainer build --fakeroot` (restricted
at many HPC sites) and a multi-GB local image build.

Keeping the install under `resources/` matters: that path is inside the working directory
and therefore already visible inside the container. An out-of-tree path still works — the
preflight appends the bind for you.

### Selecting a different location

First match wins:

```bash
snakemake --config diann_path=/opt/diann-2.5.0 ...   # 1. command line
export CONDUIT_DIANN_PATH=/opt/diann-2.5.0           # 2. environment
# 3. diann_path: in config/snakemake.yaml or your profile
```

### What gets checked, and when

| When | Check | On failure |
|---|---|---|
| At parse time | The path resolves, has one of the three shapes, and the binary is executable | Fails in seconds with download instructions. Downgraded to a warning on `--dry-run` |
| Before any search | DIA-NN's version is ≥ 2.2, and it recognises every CLI flag the workflow passes | Fails with the offending flags named. Controlled by `diann_compat_check` |
| After the main search | The report carries every column the downstream readers select by name | Fails naming the missing columns, instead of an opaque error inside an R script |

The flag check exists because DIA-NN does **not** error on an unknown flag — it prints
`WARNING: unrecognised option [--flag]` and carries on. An unchecked version mismatch
therefore changes the science silently rather than crashing. The probe runs DIA-NN with
every flag the workflow uses and no input files (a few seconds), and includes a deliberately
invalid flag as a self-test: if DIA-NN does not report *that* one, the check refuses to
claim a pass.

Set `diann_compat_check: warn` to report and continue, or `off` to skip it entirely.

The DIA-NN version, mode, resolved path, and flag findings are recorded in each run's
`manifest.json` under `diann` — the version is no longer pinned by an image tag, so runs
have to capture what actually executed.

#### What these checks do not tell you

They are **basic structural checks**, and it is worth being clear about their limits.

They confirm that a given DIA-NN *runs*, that it understands every flag Conduit passes, and
that its report carries the columns Conduit reads. They say nothing about whether it
produces the **same numbers** as the version you used last. A release can change scoring,
FDR estimation, RT modelling, or quantification and still pass every check here without a
warning — the interface is identical, the science is not.

DIA-NN releases often — twelve Linux builds since 2.0 — and we do not test every one.
"Tested" means someone ran Conduit end to end on that version and was satisfied with the
result; it is a short list, maintained by hand in `DIANN_TESTED_VERSIONS`
(`modules/_shared/diann_env.py`). An untested version that passes the checks gets a
**warning, not a blessing**: it means "nothing structural is wrong", not "this is known to
be equivalent".

Two practical consequences:

- **Pin one version for the duration of a study.** Results from different DIA-NN versions
  are not automatically comparable, and `manifest.json` records which version produced each
  run so you can tell them apart after the fact.
- **Upgrading is a change worth measuring.** Both versions can be installed side by side and
  selected per run, so the comparison is cheap:

  ```bash
  snakemake --config diann_path=resources/diann/diann-2.5.0 run_name=v250 ...
  snakemake --config diann_path=resources/diann/diann-2.6.1 run_name=v261 ...
  ```

  Then diff what you actually care about — precursor and protein-group counts, taxon calls,
  quantities — before treating the new version as a drop-in replacement.

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
│   │   ├── genomes/                  # From user-provided genome FASTAs
│   │   └── database_processing/      # Shared post-processing
│   ├── diann/                        # DIA-NN identification and quantification
│   ├── annotation/                   # Protein and taxonomic annotation
│   │   ├── uniprot/                  # UniProt-based annotation
│   │   ├── genomes/                  # Bakta + UniProt annotation for genome proteins
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
│       │   ├── ms_files/             # .raw or .mzML MS files
│       │   ├── sample_annotation.txt # Sample metadata
│       │   ├── ncbi_taxa_ids.txt     # Taxon IDs (ncbi_taxonomy_id method)
│       │   ├── proteome_ids.txt      # Proteome IDs (uniprot_proteome_id method)
│       │   ├── fastq_files/          # FASTQ files (metaphlan method)
│       │   └── genome_files/         # Genome FASTA files (genomes method)
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

> For tutorial purposes you will also need files in `experiments/example/input/database_resources` and `experiments/example/input/ms_files`. Contact baynec2 directly for these files.

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
Don't install singularity/apptainer into it — Barnacle2 has Singularity installed
system-wide (`/usr/bin/singularity`, SingularityCE 4.x), so keeping it out of the
env leaves that system Singularity on `PATH` (on every node, no module needed).

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
singularity --version       # system install at /usr/bin/singularity (SingularityCE 4.x)

# Point Snakemake's caches/temp at scratch so they don't fill your quota-limited
# $HOME. SLURM exports this environment to the per-rule jobs.
export XDG_CACHE_HOME="/ddn_scratch/${USER}/.cache"
export TMPDIR="/ddn_scratch/${USER}/tmp"
mkdir -p "$XDG_CACHE_HOME" "$TMPDIR"

# --- launch the workflow ---
snakemake \
  --profile profiles/barnacle2 \
  --configfile experiments/<exp>/config/<method>.yaml
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

conduit-ascent is released under the [MIT License](LICENSE), copyright © 2025-2026
Charlie Bayne.

Third-party material copied into this repository — and the terms of the reference data
the pipeline downloads — is inventoried in
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).

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
