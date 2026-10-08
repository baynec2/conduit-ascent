# Changelog

All notable changes to conduit-ascent are recorded here. Versions follow
[Semantic Versioning](https://semver.org/). Before 1.0.0, a minor version may change
config keys or output formats; check this file before upgrading mid-study.

## [0.1.1] - 2026-10-07

### Fixed

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
