# Provenance — `ribP_elonF_profiles_refined_manually.hmm`

77 profile HMMs for ribosomal proteins and elongation factors, used by the
`hapid` and `unipept_hapid` search-space methods to find marker genes in
genome FASTAs.

| | |
|---|---|
| Source | [Pfam](https://www.ebi.ac.uk/interpro/) (Pfam-A), via InterPro |
| License | [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/) — public-domain dedication |
| Format | HMMER3/f, built by Pfam; profiles carry their `PF#####.#` accessions |
| Selection | The accession list is the marker-gene set described by [HAPiID](https://github.com/mgtools/HAPiID) |

## What was taken from where

The **models** are Pfam's, redistributed under CC0, which places no condition on
redistribution. The **choice** of which 77 Pfam families constitute the ribosomal-protein
and elongation-factor marker set follows HAPiID's published method — a list of accessions,
not copyrightable expression, and not HAPiID's model files.

This matters because `mgtools/HAPiID` publishes no license, so nothing of theirs can be
redistributed. Rebuilding the profile set from Pfam rather than copying HAPiID's own HMM
file keeps this repository clear of that problem.

To list the accessions in the file:

```bash
grep '^ACC' resources/hapid/ribP_elonF_profiles_refined_manually.hmm
```

## Attribution

Pfam asks to be cited. See <https://www.ebi.ac.uk/interpro/about/citing/>.

The marker-gene approach is HAPiID's; credit it as well:

> Stamboulian M, Li S, Ye Y. *Using high-abundance proteins as guides for fast and
> effective peptide/protein identification from human gut metaproteomic data.*
> Microbiome 9, 80 (2021). doi:[10.1186/s40168-021-01035-8](https://doi.org/10.1186/s40168-021-01035-8)
