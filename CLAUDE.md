# conduit-ascent — Project Guide for Claude

## What This Project Is

conduit-ascent is a Snakemake workflow for DIA (Data-Independent Acquisition) metaproteomics. It builds taxonomic protein databases from raw mass spectrometry files and produces quantified, annotated protein/peptide matrices. The pipeline selects organisms present in a sample before downloading full proteomes, avoiding the prohibitive cost of searching against all of UniProt.

## Repository Structure

```
modules/
├── search_space/          # How organisms are selected (one method per subdirectory)
│   ├── unipept_peptidotyping/  # UMGAP-based tiered search (see below)
│   ├── ncbi_taxonomy/     # User-supplied NCBI taxon IDs → UniProt proteomes
│   ├── uniprot_proteome_ids/  # User-supplied UniProt proteome IDs directly
│   ├── metaphlan/         # MetaPhlAn profiling output → NCBI taxon IDs
│   ├── genomes/          # Any bacterial genome FASTAs (MAGs or reference), no selection; "MAGs" is a deprecated alias
│   └── unipept_hapid/     # Unipept database search (method: unipept_hapiid; "unipept_hapid" is a deprecated alias, dir name kept)
├── diann/                 # DIA-NN spectral library build + main search
├── annotation/            # UniProt, eggNOG-mapper, external annotations
└── build_conduit/         # Final conduit object construction (uses conduitR)
config/                    # DIA-NN and tool config files
experiments/               # Per-experiment input dirs (ms_files/*.raw|*.mzML, config YAML)
tests/                     # Integration tests
```

## Peptidotyping — Approach and Design

The `unipept_peptidotyping` method is the flagship search-space method. It identifies which organisms are present in a sample before downloading proteomes, using UMGAP-derived LCA (Lowest Common Ancestor) peptides from all of UniProt (SwissProt + TrEMBL).

### Core Idea

Every tryptic peptide in UniProt is assigned its LCA taxon by UMGAP. Peptides unique to a single family (or genus, or species/strain) serve as diagnostic markers. A tiered DIA-NN search against these diagnostic peptides identifies organisms present, producing NCBI taxonomy IDs that feed the downstream `ncbi_taxonomy` module.

### Two-Pass Tiered Search

**Pass 1 — Family-level detection:**
- Database: `effective_first_pass_database.fasta`
  - For each family: includes peptides at the finest rank (family → genus → species/strain) where ≥ `min_taxon_db_peptides` (default 10) unique peptides exist
  - Every FASTA header carries a `FAM=<family_taxid>` tag so the family is always recoverable regardless of which rank's peptides were used
- DIA-NN is run in **InfiniDIA mode** (`--pre-search --pre-filter`) with `--qvalue 1 --report-decoys` (all PSMs reported for FDR calculation)
- Family presence is determined by **picked target-decoy FDR** (Savitski et al. 2015) at the taxon level via `conduitR::calc_taxon_fdr`: it aggregates per-`(taxon, decoy)` scores, keeps only the higher-scoring of each lineage's `{target, decoy}` pair, and competes the representatives in one ranked list. A lineage is called present iff it is a target representative AND `qvalue ≤ qvalue_threshold` (default 0.05) AND `n_unique_peptides_all ≥ min_peptides` (default 2). This stops a high-abundance family's reversed decoy from outranking a true low-abundance family's target.

**Pass 2 — Species/strain resolution:**
- Database: `second_pass_database.fasta`
  - Filtered subset of `species_strain_lca_filtered_peptides.tsv` containing only taxa that are descendants of the families detected in pass 1
- Same DIA-NN InfiniDIA search
- Species/strain presence determined by the same picked target-decoy FDR (`calc_taxon_fdr`) on `OX=<taxid>` extracted from FASTA headers

**Output:** `ncbi_taxa_ids.txt` — detected species/strain NCBI taxon IDs, handed to the `ncbi_taxonomy` module which fetches UniProt proteomes.

### Key Files

| File | Purpose |
|------|---------|
| `modules/search_space/unipept_peptidotyping/unipept_peptidotyping.smk` | All Snakemake rules |
| `scripts/infer_family_presence.R` | Pass 1 picked-FDR inference (extracts `FAM=` tag) |
| `scripts/infer_species_strain_presence.R` | Pass 2 picked-FDR inference (extracts `OX=` tag) |
| `scripts/presence_lib.R` | Shared helper `apply_picked_presence_filter` (picked FDR + optional coverage filter) |
| `scripts/unipept-database/scripts/generate_umgap_tables.sh` | Builds UMGAP sequence index from UniProt |
| `config/peptidotyping_infinidia.cfg` | DIA-NN InfiniDIA config (shared by both passes) |

### Rule Dependency Graph

```
build_sequence_index  (UniProt → sequences.tsv.lz4 + taxons.tsv.lz4)
        │
generate_peptidotyping_db  (×3 ranks: family, genus, species_strain)
        │
build_effective_detection_rank_db  (→ effective_first_pass_database.fasta + rank_mapping.tsv)
        │
perform_first_pass_search  [DIA-NN InfiniDIA]
        │
infer_first_pass_presence  [calc_taxon_fdr picked FDR @ q≤0.05, extracts FAM= tag]
        │
map_first_pass_detected_taxa_to_species_strains  [taxonkit]
        │
generate_second_pass_db  [awk filter species_strain TSV]
        │
perform_second_pass_search  [DIA-NN InfiniDIA]
        │
infer_second_pass_presence  [calc_taxon_fdr picked FDR @ q≤0.05, extracts OX= tag]
        │
generate_peptidotyping_ncbi_taxa_ids  →  ncbi_taxa_ids.txt
```

### FASTA Header Formats

- **First-pass (effective_first_pass_database.fasta):**
  `umgap|{id}|{lca_il} {rank}_{name} OS={name} OX={lca_il} RK={rank} PT={parent_id} FAM={family_taxid}`
  → Parse `FAM=` for family taxon assignment

- **Second-pass (second_pass_database.fasta, from species_strain TSV):**
  `umgap|{id}|{lca_il} species_strain_{name} OS={name} OX={lca_il} RK={rank} PT={parent_id}`
  → Parse `OX=` for species/strain taxon assignment

### `calc_taxon_fdr` (conduitR)

Signature: `calc_taxon_fdr(pep, taxon, decoy, peptide = NULL, qvalue_threshold = 0.05, min_peptides = 2)`

- **Picked target-decoy FDR** (Savitski et al. 2015) in one call: aggregates PSM scores to one target row + one decoy row per taxon (`sum(-log(PEP))` + distinct-peptide counts), keeps only the higher-scoring of each lineage's `{target, decoy}` pair as the representative, discards the loser, then competes all representatives in one descending-score list
- `fdr = (cum_decoy_winners + 1) / max(cum_target_winners, 1)`; `qvalue` = cumulative min of `fdr` from the bottom up
- A lineage is **present** (`pass == TRUE`) iff it is a target representative AND `qvalue ≤ qvalue_threshold` AND `n_unique_peptides_all ≥ min_peptides`
- Mixed ranks (family/genus/species/strain) compete together — picking sinks the null to the bottom, so no per-rank stratification is needed
- Fixes the abundance bias of a running-FDR-over-both-rows scheme, where a high-abundance lineage's reversed decoy could outrank a true low-abundance lineage's target and wrongly reject it
- Returns a list with `$results` (per-`(taxon, decoy)` rows incl. `picked_winner` / `fdr` / `qvalue` / `pass`), `$detected`, `$n_targets`, `$n_decoys`, `$first_decoy_rank`, `$n_missing_pair`
- Requires unfiltered DIA-NN output (`--qvalue 1`) and decoys (`--report-decoys`); pools PSMs across all runs for experiment-level taxon calls
- The picked competition core is the internal `conduitR:::pick_taxon_fdr_compete` (operates on an aggregated score table — the unit/regression test target). In peptidotyping it is wired in via `apply_picked_presence_filter` (in `presence_lib.R`), which builds the audit table and applies the optional, **off-by-default** score-coverage filter

### Container Strategy

All container images are pinned to short-SHA tags in `config/snakemake.yaml` under the `containers:` block. Two CI workflows publish them:

- **conduitR** (`baynec2/conduitr`): published by `conduitR/.github/workflows/docker-publish.yml` on every push to `main`/`develop`. Tags emitted: `:latest` (main), `:develop` (develop), and `:<short-sha>` (immutable).
- **In-repo containers** (diann, bakta, metaphlan, eggnogmapper, umgap, fraggenescan_hmmer): published by `.github/workflows/build-container-images.yml` on push to `main`/`develop` when a `containers/*/Dockerfile` changes. Only the immutable `:<short-sha>` tag is emitted (no rolling tag).

**Pin to a SHA, not a moving tag like `:develop`/`:alpha`.** Apptainer caches images by URI in `.snakemake/singularity/`, so a moving tag won't auto-refresh once cached — you'd silently keep running an old image. A SHA pin makes the URI change explicit so the new image gets pulled. Bump intentionally when you want an upstream change.

### Pre-flight: bump container tags before starting a workflow

During active dev, container Dockerfiles change often. Before kicking off any non-trivial workflow run, verify the tags in `config/snakemake.yaml` match the latest CI-published SHA for each image:

**Tag length matters: pin to the 7-char SHA, not 8.** CI (`build-container-images.yml`) tags every image with `git rev-parse --short HEAD`, which abbreviates to **7 characters** (e.g. `e46d512`). Git's `%h` / `git log` can auto-abbreviate to 8+ chars on this repo, so always force `--abbrev=7` when reading the SHA to pin — an 8-char pin like `e46d512a` will `MANIFEST_UNKNOWN` because that tag was never pushed.

```bash
# Latest in-repo Dockerfile commit per container (7-char, matches the CI-pushed tag)
for tool in bakta diann eggnogmapper fraggenescan_hmmer metaphlan umgap; do
  printf "%-22s %s\n" "$tool" "$(git log develop --abbrev=7 --format='%h' -1 -- "containers/${tool}/Dockerfile")"
done

# Latest CI-published runs (the SHA you should pin to is the latest "success" headSha)
gh run list --workflow=build-container-images.yml --limit 5 --json conclusion,headSha,displayTitle

# Latest conduitR develop SHA (for the conduitr container)
git -C /home/nanopore-catalyst/conduitR log --oneline --abbrev=7 develop -1

# Verify a tag actually exists on Docker Hub before pinning (replace tool/sha):
curl -s "https://hub.docker.com/v2/repositories/baynec2/umgap/tags?page_size=25" | python3 -c "import sys,json;print('\n'.join(t['name'] for t in json.load(sys.stdin)['results']))"
```

Bump any tag in `config/snakemake.yaml` whose SHA differs from the latest CI-published one. (Once things stabilize and Dockerfiles aren't changing, this becomes a rare check.)

If a Dockerfile commit predates the CI workflow being added (`184ad867`, 2026-04-06) and has no SHA-tagged image on Docker Hub — bakta, eggnogmapper, metaphlan, fraggenescan_hmmer at time of writing — either leave the existing manual tag (`:alpha`, `:2.1.12`) or trigger a CI rebuild by making a no-op change to the Dockerfile.

## conduitR Package

Located at `/home/nanopore-catalyst/conduitR`. An R package providing:
- S4 `conduit` class for integrated metaproteomics data
- `calc_taxon_fdr()` — picked taxonomic target-decoy FDR (Savitski et al. 2015), used by peptidotyping
- Taxonomy utilities: `add_taxonomy_to_qf`, `get_ncbi_taxonomy`, etc.
- Quantification: `calc_relative_abundance`, `add_log_imputed_norm_assay`, etc.
- Visualization: `plot_taxa_tree`, `plot_volcano`, `plot_sunburst`, etc.

## DIA-NN Configuration Notes

- `--mass-acc 10` / `--mass-acc-ms1 4` — correct for this instrument; **do not change**
- Peptidotyping uses `--cut "" --missed-cleavages 0` because peptides are pre-digested (UMGAP tryptic digest)
- InfiniDIA flags: `--pre-search --pre-filter` (required for the broad database search strategy)

## Running the Workflow Locally

Use the per-hostname Snakemake profile under `profiles/<hostname>/`. The profile injects `--cores`, `--use-singularity`, the `/home/nanopore-catalyst/HDD` bind, and any per-machine resource-path overrides — no manual flags required:

```bash
conda activate conduit
snakemake --profile profiles/<hostname> \
  --configfile experiments/<name>/config/<method>.yaml <target>
```

`tests/run_integration_tests.sh` auto-applies the profile that matches `hostname`.

### Peptidotyping Resources

Pre-built peptidotyping resources (UMGAP sequence index, peptide TSVs, taxonkit DB) live on the HDD at `/home/nanopore-catalyst/HDD/peptidotyping_resources/`. The `nanopore-catalyst` profile already points `peptidotyping_resource_dir` and `taxonkit_db_dir` at this location, so experiment configs do **not** carry these paths. To enable a new host, copy `profiles/nanopore-catalyst/` to `profiles/<new-hostname>/` and edit the bind / resource paths to match that machine's storage layout.
