# Third-Party Notices — conduit-ascent

Third-party material **that is copied into this repository** and therefore travels with
it. Tools this pipeline merely *invokes* — DIA-NN, eggNOG-mapper, Bakta, MetaPhlAn,
TaxonKit and the rest — are not listed here: each is obtained by the user from its own
publisher, runs as a separate program inside its own container, and is bound by its own
terms rather than by anything in this file. See the "Dependencies" section of
[`README.md`](README.md) for what the pipeline expects and downloads at run time.

## Code

### unipept-database — MIT

`modules/search_space/unipept_peptidotyping/scripts/unipept-database/`

An unmodified partial copy of <https://github.com/unipept/unipept-database> at commit
`7f4e5951e25c`, copyright (c) 2023 Universiteit Gent. MIT permits the copy and requires
its notice to travel with it.

- License: [`.../unipept-database/LICENSE`](modules/search_space/unipept_peptidotyping/scripts/unipept-database/LICENSE)
- Provenance and verification: [`.../unipept-database/PROVENANCE.md`](modules/search_space/unipept_peptidotyping/scripts/unipept-database/PROVENANCE.md)

### coverAllSpectra_greedy.py — derived from HAPiID, no license declared

`modules/search_space/_shared/scripts/coverAllSpectra_greedy.py`

Adapted from <https://github.com/mgtools/HAPiID> (Ye lab, Indiana University), which
publishes **no license file** — so default copyright applies upstream and this file
cannot be relied on as redistributable. Recorded here rather than left implicit.

The underlying greedy set-cover algorithm is textbook and not protectable; what remains
at issue is the specific expression of two helper functions and their naming. Resolving
it needs either written permission from the upstream authors or an independent
reimplementation of those helpers. **Open — see the licensing brief, item 2.**

## Data

### Pfam profile HMMs — CC0 1.0

`resources/hapid/ribP_elonF_profiles_refined_manually.hmm`

77 Pfam-A profile HMMs, public-domain dedication. The accession set follows HAPiID's
marker-gene method; the models themselves are Pfam's, not HAPiID's.

- Provenance: [`resources/hapid/PROVENANCE.md`](resources/hapid/PROVENANCE.md)

### Test fixtures and example inputs

`tests/fixtures/peptidotyping_subset/*.tsv` are small slices of the UMGAP LCA peptide
tables: UniProtKB peptide sequences (CC BY 4.0) carrying GO (CC BY 4.0), InterPro (CC0)
and EC cross-references, keyed to NCBI Taxonomy identifiers (public domain).
`tests/fixtures/hapid_subset/` holds the same Pfam HMM file described above.

`experiments/integration_test/input/` holds NCBI-derived taxonomy and MetaPhlAn clade
tables (public domain / MIT) alongside lab-generated `.mzML` and `.fastq.gz` inputs.

All of it is attribution-only at most, and attribution for the redistributed portion is
given by this file.

## Reference data the pipeline downloads

Not redistributed here, but their content ends up in pipeline output, so their terms
follow the results:

| Resource | Terms |
|---|---|
| UniProtKB (SwissProt + TrEMBL) | CC BY 4.0 — attribution required |
| NCBI Taxonomy | US Government public domain |
| Gene Ontology | CC BY 4.0 — attribution must carry creator, copyright, license and disclaimer |
| Pfam / InterPro | CC0 1.0 |
| ExPASy ENZYME | CC BY 4.0, © SIB Swiss Institute of Bioinformatics |
| eggNOG 5.0 database | CC BY 3.0 (note: the eggNOG-mapper *tool* is AGPL-3.0) |
| MetaPhlAn / ChocoPhlAn | MIT |
| Bakta database | Aggregate; inherits the terms of UniProt, RefSeq and its other sources |
| KEGG | Proprietary. Free for academic use; commercial use needs a Pathway Solutions license. Read live, not redistributed |
| CAZy / dbCAN | © 1998–2026 AFMB–CNRS–AMU–INRAE. Citation guidance only. Read live, not redistributed |

## DIA-NN

DIA-NN is required, proprietary, and **not distributed by this project**. Each user
downloads it from the vendor and points the workflow at their own copy via `diann_path`,
accepting the vendor's terms directly. See "Obtaining DIA-NN" in [`README.md`](README.md).
