# Changelog

All notable changes to conduit-ascent are recorded here. Versions follow
[Semantic Versioning](https://semver.org/). Before 1.0.0, a minor version may change
config keys or output formats; check this file before upgrading mid-study.

## [0.1.1] - 2026-11-17

### Added

- Optional MS run subsets for the DIA-NN steps that learn from the runs rather than
  quantify them, to bound memory on large cohorts. Two new config keys, both `null`
  (every run) by default, so existing configs behave as before:
  `diann_empirical_lib_subset` builds the empirical library from a subset (stages 2–3
  still search and quantify every run), and `search_space_run_subset` runs the
  organism-selection passes of `unipept_peptidotyping`, `unipept_hapiid`, `hapiid` and
  `genome_peptidotyping` on a subset. Each takes `{fraction, stratify_by, seed}` or
  `{samples: [...]}`; the selection is recorded in `ms_subsets/<name>.selection.tsv`
  and `manifest.json`. An organism present only in left-out runs never enters the
  database, so prefer a generous fraction. Subsets are staged as hard links, because
  DIA-NN 2.5.0 cannot read a Thermo `.raw` through a symlink.

### Fixed

- The staged DIA-NN combine step no longer re-searches every run. The workflow
  mispredicted DIA-NN's `.quant` file names, and a mismatch is not an error: DIA-NN
  silently searched each run again instead of reusing the per-run results. Combine
  now passes `--quant-ori-names`. Results are unchanged; the per-run stage simply
  stops being wasted work. Completed combine steps rerun once after upgrading.

- The `conduitr` container is re-pinned to conduitR 0.1.0. v0.1.0 ran a conduitR
  build from July, which predates the `get_ncbi_taxonomy()` fix: for a proteome
  registered under a strain taxid, the strain name was written into `species`
  (e.g. "Akkermansia muciniphila ATCC BAA-835"). Two strains of one species then
  had different `species` values, so peptides shared between them fell to genus and
  dropped out of species-level results. Taxonomy tables now carry the species in
  `species` and the strain designation in a separate `strain` column. Affects the
  `uniprot_proteome_id`, `ncbi_taxonomy_id`, `metaphlan`, `unipept_peptidotyping` and
  `unipept_hapiid` methods, which all build their taxonomy through it, plus `genomes`
  runs that append extra proteomes or taxa. Re-run affected v0.1.0 results
  if you rely on species-level aggregation.

### Changed

- The release checklist in `CLAUDE.md` now re-pins `conduitr` to the latest conduitR
  release, and documents how to read an image's digest.

## [0.1.0] - 2026-10-07

First tagged release.

### Search-space methods

- `ncbi_taxonomy_id` and `uniprot_proteome_id`: build the search space from user-supplied
  NCBI taxon IDs or UniProt proteome IDs
- `unipept_peptidotyping`: two-pass tiered DIA-NN search against UMGAP-derived
  diagnostic peptides, with presence calls from picked target-decoy FDR at the taxon level
- `unipept_hapiid` and `hapiid`: HAPiID-style selection of the smallest genome or
  set that covers the annotated spectra
- `genome_peptidotyping`: peptidotyping against peptides digested from your own genomes
- `metaphlan`: search space from MetaPhlAn profiling of shotgun metagenomes
- `genomes`: Bakta-annotated genome FASTAs, supplied by hand or downloaded from an
  MGnify catalog

### Workflow

- DIA-NN 2.2 or newer, downloaded by the user; version and flag compatibility are checked
  before any search
- Taxonomic and functional annotation (UniProt, eggNOG-mapper, GO, KEGG, Pfam, CAZy)
- Final `conduit` R object for conduitR and conduit-summit
- Every run writes `manifest.json` with the workflow version, git commit, resolved
  config, input hashes and container tags. Untagged commits record their version as
  `0.1.0+<commit>`, so development runs cannot pass for a release
- Containers pinned to immutable image tags
- Snakemake profiles for a local machine and the Barnacle SLURM cluster
- Released under the MIT License

[0.1.1]: https://github.com/baynec2/conduit-ascent/releases/tag/v0.1.1
[0.1.0]: https://github.com/baynec2/conduit-ascent/releases/tag/v0.1.0
