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

### coverAllSpectra_greedy.py — our code, HAPiID's method

`modules/search_space/_shared/scripts/coverAllSpectra_greedy.py`

**No third-party code.** The greedy genome-selection step implements the method
described in HAPiID, and the paper is credited in the script and below, but the
implementation is ours.

It did not start that way. The file began as an adaptation of
<https://github.com/mgtools/HAPiID>, which declares **no license at all** — so default
copyright applies upstream and nothing from it can be redistributed. Two helper
functions still carried upstream expression. They were replaced with an independent
implementation written against the algorithm and this module's tests, verified
equivalent by fuzzing 50,000 randomised inputs against the previous behaviour with zero
mismatches. The greedy set-cover algorithm itself is textbook and not protectable.

Credit for the method belongs upstream:

> Stamboulian M, Li S, Ye Y. *Using high-abundance proteins as guides for fast and
> effective peptide/protein identification from human gut metaproteomic data.*
> Microbiome 9, 80 (2021). doi:[10.1186/s40168-021-01035-8](https://doi.org/10.1186/s40168-021-01035-8)

Asking the authors to add a license to their repository remains worth doing — it would
settle the question for every other group building on HAPiID, several of whom will hit
the same wall — but this repository no longer depends on the answer.

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
