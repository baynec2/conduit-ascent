# Changelog

All notable changes to conduit-ascent are recorded here. Versions follow
[Semantic Versioning](https://semver.org/). Before 1.0.0, a minor version may change
config keys or output formats; check this file before upgrading mid-study.

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

[0.1.0]: https://github.com/baynec2/conduit-ascent/releases/tag/v0.1.0
