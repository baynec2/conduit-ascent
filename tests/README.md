# Tests

Three-tier test architecture. Pick the tier matching what your change touches.

| Tier | What | Where | When | How to run |
|---|---|---|---|---|
| **1. Unit** | Pure-function script logic on synthetic data. No DBs, no DIA-NN, no containers. | `tests/scripts/{testthat,python}/`, `tests/scripts/test_awk_pipelines.sh`, `tests/scripts/test_parse_eggnogmapper_annotations.R` | CI, every push | `Rscript -e "testthat::test_dir('tests/scripts/testthat')"` · `pytest tests/scripts/python` · `bash tests/scripts/test_awk_pipelines.sh` |
| **2. Smoke** | End-to-end pipeline against subset DBs. Heavy enough to catch real DIA-NN/snakemake integration bugs, light enough to iterate. | `tests/run_smoke_tests.sh`, `tests/configs/smoke_*.yaml`, `tests/fixtures/` | Local-only, run pre-push | `bash tests/run_smoke_tests.sh [method]` |
| **3. Full integration** | End-to-end pipeline against production DBs (~170 GB UMGAP + bakta + eggNOG). | `tests/run_integration_tests.sh`, `tests/configs/integration_test_*.yaml` | Local, before releases / on demand | `bash tests/run_integration_tests.sh [method]` |

## Tier 1 — Unit tests

Test the **glue logic in snakemake scripts** — regex extraction, dataframe joins, threshold filters, format-parsing edge cases. Not the conduitR package internals (those have their own testthat suite at `conduitR/tests/testthat/`).

Refactor pattern when adding a new unit test for an existing script:

1. Factor the glue logic out of the snakemake-bound script into a pure function in a sibling `*_lib.R` / extracted top-level Python function.
2. Have the snakemake script source/import that pure function and call it.
3. Write the test against the pure function.

Snakemake R scripts source siblings via `snakemake@source("./file_lib.R")`. Python scripts gate their snakemake-bound `main()` on `if __name__ == "__main__" or "snakemake" in globals():` so plain `import` runs no I/O.

### Currently covered

| Test file | Function under test | Located in |
|---|---|---|
| `testthat/test-presence_lib.R` | `normalize_param`, `extract_family_psms`, `extract_species_strain_psms` | `modules/search_space/unipept_peptidotyping/scripts/presence_lib.R` |
| `testthat/test-call_ncbi_taxa_ids_lib.R` | `parse_metaphlan_profiles` | `modules/search_space/metaphlan/scripts/call_ncbi_taxa_ids_lib.R` |
| `testthat/test-get_uniprot_proteome_ids_lib.R` | `parse_organism_ids` | `modules/search_space/ncbi_taxonomy/scripts/get_uniprot_proteome_ids_lib.R` |
| `test_parse_eggnogmapper_annotations.R` | inline regression cases | `modules/annotation/eggnogmapper/scripts/parse_eggnogmapper_annotations.R` |
| `python/test_parse_genome_taxonomy.py` | `build_genome_taxonomy`, `assign_organism_ids`, `filter_to_selected_genomes`, `fill_taxonomy_defaults` | `modules/search_space/genomes/scripts/parse_genome_taxonomy.py` |
| `python/test_build_taxon_spectrum_mapping.py` | `build_taxon_spectrum_mapping` | `modules/search_space/unipept_hapid/scripts/build_taxon_spectrum_mapping.py` |
| `python/test_cover_all_spectra_greedy.py` | `greedy_cover` | `modules/search_space/_shared/scripts/coverAllSpectra_greedy.py` |
| `test_awk_pipelines.sh` | FASTA header generation + rank-priority awk pipelines from `unipept_peptidotyping.smk` | `modules/search_space/unipept_peptidotyping/unipept_peptidotyping.smk` |

### Gaps

- `modules/search_space/genomes/scripts/genome_uniprot_headers.py` — needs bakta `.tsv` / `.faa` fixtures; deferred.
- Other R/python scripts in `modules/` don't yet have a `*_lib.R` factoring.
- AWK tests duplicate the pipeline code from the `.smk`; if the `.smk` changes, the test must be updated in lockstep (no shared source).

## Tier 2 — Smoke tests

Smoke tests **run the real DIA-NN + downstream pipeline** but against subset fixtures so the cycle is minutes, not hours. They exist primarily to catch:

- Format mismatches at module seams that unit tests can't see.
- DIA-NN config regressions (precursor m/z windows, FDR settings, etc.).
- Container / Singularity / apptainer plumbing breakage.

### Current smoke coverage

| Method | Config | Fixture | Detects |
|---|---|---|---|
| `unipept_peptidotyping` | `tests/configs/smoke_unipept_peptidotyping.yaml` | `tests/fixtures/peptidotyping_subset/` (10 pool species + ancestors, ~17 MB) | 4 of 10 pool families against the pool test mzML |

The fixture peptide TSVs are LFS-tracked. They can be regenerated from the full UMGAP index via `tests/scripts/build_subset_peptidotyping_db.sh`.

### Why smoke is local-only

The peptidotyping subset fixture is 17 MB committed via LFS; the test mzML is 1.1 GB via LFS. CI on a GitHub-hosted runner would need to LFS-pull both per run (slow + LFS bandwidth) and pull multiple ~2 GB singularity images. Local smoke is the practical path; CI stays on dry-runs + Tier 1 units.

### Adding a smoke method

1. Create `tests/configs/smoke_<method>.yaml` modelled after `smoke_unipept_peptidotyping.yaml`. Override `peptidotyping_resource_dir` / `hapiid_hmm_profiles` etc. to fixture paths if the method needs them.
2. Add a `run_<method>` function to `tests/run_smoke_tests.sh` and wire it into the dispatch `case`.
3. If the method needs a fixture, build it under `tests/fixtures/<method>_subset/` and LFS-track it via `.gitattributes`.

## Tier 3 — Full integration tests

The original integration suite at `tests/run_integration_tests.sh` against production DBs. Use this when:

- Validating a release candidate end-to-end.
- Diagnosing a problem that only shows up at full DB scale (e.g., FDR statistics).
- Updating an immutable container tag.

Most methods take several minutes each; the peptidotyping family takes longer because of DIA-NN against the full UMGAP index. eggNOG annotation runtime varies dramatically by input dataset (the `--iterate` flag does multiple sensitivity passes if hit count is low).

## Adding a new test — which tier?

| Change | Tier |
|---|---|
| New regex, dataframe join, threshold logic, format parser | **Unit** (extract pure function + testthat / pytest) |
| New snakemake rule / DAG edge / shell pipeline | **Smoke** if it can run against the subset fixture; **integration** otherwise |
| New search-space method | All three: unit tests for the new scripts; smoke config with subset fixture; integration config with full DB |
| Container image bump | At minimum a smoke run; integration if the rule's logic depends on container content |
| Config schema change | Smoke + the relevant integration variants |

## Known caveats

- **CI runs only Tier 1 + dry-runs**, not smoke. Local smoke is the contract for catching pipeline regressions before push.
- **The pool test mzML (1.1 GB) replaced the old 134 MB E. coli mzML.** Integration tests that previously detected only E. coli will now match against pool species too. If a method's protein detection count differs significantly between methods, that's expected — they have different upstream taxon selection.
- **eggNOG annotation runtime is dataset-dependent.** Diamond's `--iterate` mode can run 10× longer on Swiss-Prot vs TrEMBL inputs even when protein counts match. This is annotation-step behavior, not a workflow regression. The Tier 1 + smoke layers don't depend on eggNOG, so this doesn't affect day-to-day iteration.
- **`conduitR::get_better_proteome_ids()`** crashes with an opaque `dplyr::left_join` error on empty input (when MetaPhlAn detects no taxa above threshold). Would benefit from a guard upstream. Not blocking current tests because we ensure MetaPhlAn detects at least one species.
