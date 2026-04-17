# HAPiID Module Implementation Plan

## Context

This plan implements `hapid` as a new `search_space_method` in conduit-ascent. HAPiID is a marker-gene-based organism profiling approach that identifies which genomes are present in a metaproteomics sample, then builds a targeted protein database from those genomes for a final DIA-NN search.

**Reference implementation:** https://github.com/mgtools/HAPiID — used as the reference for all pipeline logic (FGS params, HMMER E-value cutoffs, CD-HIT settings, greedy selection algorithm). The HAPiID paper (Microbiome, 2021) describes the algorithm; the paper originally uses MSGF+ as the search engine — this plan adapts it for DIA-NN to support DIA data.

**Template within conduit-ascent**: The MAG module (`modules/search_space/MAGs/MAGs.smk`) — genome-FASTA-in, Bakta-annotated, UniProt-style-database-out. HAPiID adds a marker-gene profiling and greedy selection stage *before* Bakta so that Bakta only runs on the selected subset of genomes.

**Key design choices:**
- Two input modes: user-supplied genome FASTAs (custom) or precompiled marker gene database
- User provides a taxonomy file (`hapid_taxonomy.txt`) with full lineage per contig (domain → species) — no NCBI ID lookup required; enables downstream taxonomic aggregation
- Only ONE new container: `fraggenescan_hmmer` (FGS + HMMER + CD-HIT); all other rules reuse existing `config["containers"]["diann/bakta/conduitr"]`
- `mag_annotation` module reused verbatim — no new annotation module needed
- Non-selected genomes intentionally excluded from final database (conduit requires full proteome coverage per organism for taxonomic aggregation)

---

## Two Modes

### Mode A — Custom database (user provides genome FASTAs)

User provides:
```
experiments/{experiment}/input/hapid_genome_files/
├── hapid_taxonomy.txt     # contig_id  domain  kingdom  phylum  class  order  family  genus  species
├── genome_A.fasta         # contigs for genome A
└── genome_B.fna           # contigs for genome B
```
`hapid_taxonomy.txt` assigns a full taxonomy row to every contig (TSV, header row required). Each FASTA file = one genome/MAG/bin. The genome ID is the FASTA basename (no extension).

### Mode B — Precompiled database

User sets `hapid_use_precompiled_db: true` and `hapid_precompiled_db: "resources/hapid/precompiled_marker_genes.fasta"` in config. Steps 1–4 (FGS, HMMER, CD-HIT, build marker gene fasta) are skipped entirely. The precompiled FASTA must already have `{genome_id}|{protein_id}` headers consistent with the taxonomy file provided.

The taxonomy file is still required in mode B (needed for Bakta annotation and conduit taxonomy.txt).

---

## Full Pipeline Shape

```
[Mode A only]
All input genomes (*.fa / *.fna / *.fasta)
    │
    ├─ check_hapid_fastas
    ├─ predict_orfs_with_fraggenescan       ─ scatter by {genome}
    ├─ identify_marker_genes_with_hmmer     ─ scatter by {genome}
    ├─ build_hapid_marker_gene_fasta          aggregate HMMER hits → all_marker_genes.fasta
    ├─ deduplicate_marker_genes_with_cdhit    CD-HIT -c 1.0 → marker_gene_db.fasta + .clstr
    │                                         also write: protein2genome_dic.json
[Mode B: inject precompiled_db here as marker_gene_db.fasta]
    │
    ├─ create_hapid_profiling_spectral_library   (DIA-NN predict)
    ├─ perform_hapid_profiling_search            (DIA-NN --dir → profiling_report.parquet)
    │
    ├─ build_genome_spectrum_mapping       parquet + .clstr → genome2spectrum_dic.json
    ├─ [CHECKPOINT] run_greedy_genome_selection → hapid_greedy_selection.tsv
    │               ↓ dynamic genome list
    ├─ press_hapid_hmm_profiles (prereq, idempotent)
    ├─ annotate_selected_hapid_genomes_with_bakta   ─ scatter by {genome} (selected only)
    ├─ create_hapid_uniprot_style_database     full proteomes of selected genomes only
    ├─ parse_hapid_taxonomy                    contig taxonomy file → taxonomy.txt
    ├─ get_hapid_annotations  →  bakta/mag_annotations.txt  (same path mag_annotation expects)
    └─ append_hapid_additional_organisms_or_proteomes
            ↓
    database.fasta + taxonomy.txt
            ↓
    database_processing → diann → mag_annotation (reused) → external_annotation → build_conduit
```

---

## Detailed Rule Spec

### Config additions (`config/snakemake.yaml`)

```yaml
# --- HAPiID ---
hapid_use_precompiled_db: false
hapid_precompiled_db: ""                   # path if hapid_use_precompiled_db: true
hapid_hmm_profiles: "resources/hapid/ribP_elonF_profiles_refined_manually.hmm"
hapid_percent_spectra: 80                  # paper default: genomes covering 80% of profiling spectra; user-adjustable
hapid_fgs_hmmer_container: "docker://baynec2/fraggenescan_hmmer:alpha"
# All other tool containers reuse existing config["containers"] keys:
#   diann    → config["containers"]["diann"]
#   bakta    → config["containers"]["bakta"]
#   conduitr → config["containers"]["conduitr"]
```

### Top-level helpers (`hapid.smk`)

```python
import os, glob, pandas as pd

EXPERIMENT_DIR = config["experiment_dir"]
RUN_DIR        = config["run_dir"]
HAPID_DIR      = os.path.join(EXPERIMENT_DIR, "input/hapid_genome_files")
BAKTA_DIR      = config["bakta_db_dir"]
BAKTA_OUT_ROOT = os.path.join(RUN_DIR, "database_resources/bakta")
DB_OUT_ROOT    = os.path.join(RUN_DIR, "database_resources")
bakta_db_final = f"{BAKTA_DIR}-{config['bakta_db_type']}"

USE_PRECOMPILED = config.get("hapid_use_precompiled_db", False)

def get_all_hapid_genomes():
    genomes = []
    for ext in ("fa", "fna", "fasta"):
        for f in glob.glob(os.path.join(HAPID_DIR, f"*.{ext}")):
            genomes.append(os.path.splitext(os.path.basename(f))[0])
    return sorted(set(genomes))

def hapid_fasta_path(wildcards):
    for ext in ("fa", "fna", "fasta"):
        p = os.path.join(HAPID_DIR, f"{wildcards.genome}.{ext}")
        if os.path.exists(p): return p
    return os.path.join(HAPID_DIR, f"{wildcards.genome}.fa")

def get_selected_genomes(wildcards):
    chk = checkpoints.run_greedy_genome_selection.get(**wildcards).output[0]
    df = pd.read_csv(chk, sep="\t")
    pct = config.get("hapid_percent_spectra", 80)
    # Take all genomes up to and including the one that pushes cumulative coverage past the threshold
    cutoff = df[df["cumulative_pct"] >= pct].index[0] + 1
    return df["genome"].tolist()[:cutoff]

def marker_gene_db_input(wildcards):
    if USE_PRECOMPILED:
        return config["hapid_precompiled_db"]
    return os.path.join(DB_OUT_ROOT, "hapid/marker_gene_db.fasta")
```

---

### Stage 1 — ORF prediction + marker gene identification (Mode A only)

**`check_hapid_fastas`** — mirrors MAGs.smk:42–67; validates HAPID_DIR is non-empty; output: `.fastas_checked`.

**`predict_orfs_with_fraggenescan`** — scatter over `get_all_hapid_genomes()`:
- Input: genome FA, `.fastas_checked`
- Command: `FragGeneScan -s {genome_fa} -o {output_prefix} -w 0 -t complete`
- Output: `hapid/fgs/{genome}.faa`
- Container: `config["hapid_fgs_hmmer_container"]`

**`press_hapid_hmm_profiles`** — idempotent prerequisite:
- Command: `hmmpress {config[hapid_hmm_profiles]}`
- Output: `touch(config["hapid_hmm_profiles"] + ".pressed")`
- Container: `config["hapid_fgs_hmmer_container"]`

**`identify_marker_genes_with_hmmer`** — scatter over all genomes:
- Input: `{genome}.faa`, pressed HMM profiles
- Command: `hmmscan -E 1e-10 --tblout {output} {hmm_profiles} {faa}`
- Output: `hapid/hmmer/{genome}_hmmer.txt`
- Container: `config["hapid_fgs_hmmer_container"]`

**`build_hapid_marker_gene_fasta`** — aggregate all HMMER hits:
- Input: all `{genome}_hmmer.txt` + all `{genome}.faa`
- Script: `scripts/build_marker_gene_fasta.py`
  - Parse HMMSCAN tabular output: collect protein IDs with hits (E-value ≤ 1e-10)
  - Pull those sequences from `.faa` files
  - Write headers as `>{genome_id}|{protein_id}` (pipe-delimited — critical for mapping step)
  - Remove premature stop codons (`*`) from sequences
- Output: `hapid/all_marker_genes.fasta`

**`deduplicate_marker_genes_with_cdhit`** — explicit CD-HIT step:
- Input: `hapid/all_marker_genes.fasta`
- Command: `cd-hit -i {input} -o {output.fasta} -c 1.0 -n 5 -T {threads} -d 0 -M 30000`
- Output: `hapid/marker_gene_db.fasta`, `hapid/marker_gene_db.fasta.clstr`
- Container: `config["hapid_fgs_hmmer_container"]` (cd-hit included in same image)
- Also runs script `scripts/build_protein_genome_dict.py`:
  - Reads `marker_gene_db.fasta` headers → parse `{genome_id}|{protein_id}`
  - Outputs `hapid/protein2genome_dic.json` → used by genome spectrum mapping

---

### Stage 2 — DIA-NN profiling search

**`create_hapid_profiling_spectral_library`**:
- Input: `marker_gene_db_input(wildcards)` (mode-aware), `config/peptidotyping_firstpass_diann_spectral_library.cfg`
- Command: `diann --cfg {cfg} --fasta {fasta} --threads {threads} --out-lib {output_prefix}`
- Output: `hapid/marker_gene.predicted.speclib`
- Container: `config["containers"]["diann"]`

**`perform_hapid_profiling_search`**:
- Input: `raw_files/`, spectral library, marker_gene_db.fasta, `config/peptidotyping_firstpass_diann.cfg`
- Command: `diann --cfg {cfg} --fasta {fasta} --out {output_prefix} --dir {raw_dir} --lib {speclib} --threads {threads}`
- Output: `hapid/profiling_report.parquet`
- Container: `config["containers"]["diann"]`

---

### Stage 3 — Greedy selection (checkpoint)

**`build_genome_spectrum_mapping`** — `scripts/build_genome_spectrum_mapping.py`:
- Read `profiling_report.parquet` → `Protein.Ids` column
- Parse `{genome_id}|{protein_id}` headers to map spectra → genome
- **Expand via CD-HIT cluster file** (`marker_gene_db.fasta.clstr`): for each representative protein, credit all cluster members' genomes with the same spectra (critical — same logic as HAPiID's `get_genome2spectrumMappintgsMatrix.py`)
- Output: `hapid/genome2spectrum_dic.json` (genome → list of spectrum IDs)

**`run_greedy_genome_selection`** — declared as Snakemake `checkpoint`:
- Script: `scripts/coverAllSpectra_greedy.py` (copied verbatim from HAPiID repo)
- Input: `hapid/genome2spectrum_dic.json`
- Output: `hapid/hapid_greedy_selection.tsv` (columns: genome, nSpectraCovered, cumulative_pct)
- Container: `config["containers"]["conduitr"]` (has Python + pandas)
- Outputs ALL genomes in ranked order with cumulative coverage %; `hapid_percent_spectra` threshold applied in `get_selected_genomes()`

---

### Stage 4 — Bakta annotation on selected genomes

**`download_bakta_resources`** — identical to MAGs.smk:72–92 (duplicate; no `use rule` cross-module); uses `config["containers"]["bakta"]`.

**`annotate_selected_hapid_genomes_with_bakta`** — mirrors MAGs.smk:97–144:
- Input lambda: `get_selected_genomes(wildcards)` → genomes covering `hapid_percent_spectra`% of profiling spectra
- Wildcard: `{genome}` (from selected list)
- Output: `directory(BAKTA_OUT_ROOT/{genome})`
- Shell: identical to `annotate_mags_with_bakta` (TMPDIR dance, same bakta CLI)
- Container: `config["containers"]["bakta"]`

**`create_hapid_uniprot_style_database`**:
- Input: selected Bakta dirs (lambda), user `hapid_taxonomy.txt` (for species names)
- Script: `scripts/hapid_uniprot_headers.py` — adapted from `MAGs/scripts/MAG_uniprot_headers.py`:
  - Replace `mag` → `genome`, `MAG_metadata` → `hapid_taxonomy`
  - Reads species name from `hapid_taxonomy.txt` keyed by contig_id
  - Sets `OX={genome_basename}` (no NCBI ID lookup)
- Output: `hapid_database.fasta` (UniProt-style `tr|{locus}|{tag} ... OS={species} OX={genome_id} GN={gene}`), `go_annotations.txt`, `kegg_annotations.txt`

---

### Stage 5 — Taxonomy

**Taxonomy input format** (`hapid_taxonomy.txt` in `input/hapid_genome_files/`):
```
contig_id  domain  kingdom  phylum  class  order  family  genus  species
contig_1   Bacteria  <NA>  Proteobacteria  Gammaproteobacteria  ...  Escherichia  Escherichia coli
```
- Header row required; TSV-delimited
- `contig_id` must match sequence IDs in the genome FASTA files
- Taxonomy derives from whatever classifier was used on the input genomes (GTDB-Tk, CAT, etc.)

**`parse_hapid_taxonomy`** — new rule:
- Script: `scripts/parse_hapid_taxonomy.py` (new Python script; no R/conduitR needed)
- Logic:
  1. Read user `hapid_taxonomy.txt`
  2. For each genome (FASTA basename), scan sequence headers from that FASTA file to collect its contig IDs
  3. Take the mode species/taxonomy across all contigs in the genome (handles mixed bins gracefully)
  4. Assign `organism_id = {genome_basename}` — matches `OX=` in database headers
  5. Construct `lineage` as semicolon-joined: `{domain};{kingdom};{phylum};{class};{order};{family};{genus};{species}`
  6. Write conduit-format output: `organism_id\tspecies\tlineage\torganism_type\tdownload_info`
     (`organism_type = "hapid"`, `download_info = "user_provided"`)
- Output: `DB_OUT_ROOT/hapid_taxonomy.txt`
- Container: `config["containers"]["conduitr"]`

**`get_hapid_annotations`** — reuses `MAGs/scripts/get_mag_annotations.R` verbatim:
- Input lambda: selected Bakta dirs
- Output: `database_resources/bakta/mag_annotations.txt` (exact path `mag_annotation` module expects)
- Container: `config["containers"]["conduitr"]`

**`append_hapid_additional_organisms_or_proteomes`** — reuses `MAGs/scripts/append_additional_organisms_or_proteomes.R` verbatim:
- Inputs: `hapid_database.fasta`, `hapid_taxonomy.txt`
- Output: `database.fasta`, `taxonomy.txt`
- Container: `config["containers"]["conduitr"]`

---

## Files to Create

```
modules/search_space/hapid/
├── hapid.smk
└── scripts/
    ├── build_marker_gene_fasta.py          # HMMER hits → all_marker_genes.fasta with genome|protein headers
    ├── build_protein_genome_dict.py        # marker_gene_db.fasta headers → protein2genome_dic.json
    ├── build_genome_spectrum_mapping.py    # parquet + .clstr + protein2genome → genome2spectrum_dic.json
    ├── coverAllSpectra_greedy.py           # Copied verbatim from https://github.com/mgtools/HAPiID
    ├── parse_hapid_taxonomy.py             # User contig taxonomy file → conduit taxonomy.txt format
    └── hapid_uniprot_headers.py            # Adapted from MAGs/scripts/MAG_uniprot_headers.py

containers/fraggenescan_hmmer/
└── Dockerfile                             # micromamba: fraggenescan + hmmer + cd-hit
```

## Files to Modify

| File | Change |
|---|---|
| `Snakefile` | Add `"hapid"` to `ALLOWED_METHODS`; add `module hapid` declaration; add `hapid` activation block |
| `config/snakemake.yaml` | Add HAPiID config keys (see above) |
| `.github/workflows/dry-run-tests.yml` | Add `fraggenescan_hmmer` container build; add `hapid` to DAG loop |
| `tests/configs/integration_test_hapid.yaml` | New test config (new file) |

---

## Snakefile Activation Block

```python
if config["search_space_method"] == "hapid":
    use rule * from hapid
    use rule * from database_processing
    use rule * from diann
    use rule * from mag_annotation   # reused verbatim — no new annotation module needed
    use rule * from external_annotation
    use rule build_conduit from build_conduit
```

---

## Container (`containers/fraggenescan_hmmer/Dockerfile`)

```dockerfile
FROM mambaorg/micromamba:1.5.8
USER root
RUN mkdir -p /work && chown -R $MAMBA_USER:$MAMBA_USER /work
USER $MAMBA_USER
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    fraggenescan hmmer cd-hit python biopython pandas \
    && micromamba clean --all --yes
WORKDIR /work
CMD ["/bin/bash"]
```

---

## Critical Design Notes

1. **CD-HIT cluster expansion** in `build_genome_spectrum_mapping.py` is required for correct genome ranking. Without it, genomes sharing identical marker gene sequences will be undercounted.

2. **`{genome_id}|{protein_id}` header format** must be consistent across: `build_marker_gene_fasta.py` output, `deduplicate_marker_genes_with_cdhit` (CD-HIT preserves headers), and `build_genome_spectrum_mapping.py` parsing.

3. **`organism_id` = genome filename basename** (without extension). Used as the `OX=` value in UniProt-style headers and as the key in `taxonomy.txt`.

4. **Contig-to-genome assignment** for taxonomy: parse sequence IDs from each genome FASTA to associate contigs → genome. Contig IDs in the user taxonomy file must match sequence IDs in the FASTA files.

5. **Precompiled DB mode**: when `hapid_use_precompiled_db: true`, `check_hapid_fastas`, FGS, HMMER, `build_hapid_marker_gene_fasta`, and `deduplicate_marker_genes_with_cdhit` rules are all skipped. `marker_gene_db_input()` returns the precompiled path. No `.clstr` available — `build_genome_spectrum_mapping.py` must skip cluster expansion gracefully.

6. **Container reuse**: Only `hapid_fgs_hmmer_container` is new. All other rules use `config["containers"]["diann"]`, `config["containers"]["bakta"]`, or `config["containers"]["conduitr"]`.

7. **`mag_annotation` reuse**: works without modification because `get_hapid_annotations` outputs to `database_resources/bakta/mag_annotations.txt` — the exact path the module expects.

---

## Verification

1. **Dry-run**: `snakemake -n --configfile tests/configs/integration_test_hapid.yaml` — DAG must build cleanly through checkpoint.

2. **Integration test** with 3 small genome FASTAs + taxonomy file; confirm these intermediates appear:
   - `hapid/fgs/{genome}.faa`
   - `hapid/hmmer/{genome}_hmmer.txt`
   - `hapid/all_marker_genes.fasta`
   - `hapid/marker_gene_db.fasta` + `.clstr`
   - `hapid/protein2genome_dic.json`
   - `hapid/profiling_report.parquet`
   - `hapid/genome2spectrum_dic.json`
   - `hapid/hapid_greedy_selection.tsv`
   - `hapid_database.fasta`
   - `database.fasta` + `taxonomy.txt`

3. **Taxonomy format**: `taxonomy.txt` columns: `organism_id\tspecies\tlineage\torganism_type\tdownload_info`. `organism_id` must match `OX=` values in `database.fasta` headers.

4. **Precompiled mode**: set `hapid_use_precompiled_db: true`; confirm FGS/HMMER rules absent from DAG.

5. **Existing methods unaffected**: re-run dry-runs for all five existing methods after Snakefile changes.
