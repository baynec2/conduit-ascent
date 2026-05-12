#!/usr/bin/env python3
"""
build_taxon_spectrum_mapping.py

Map DIA-NN first-pass results (parquet) to species/strain taxa for greedy
HAPiID-style set-cover selection.

Each peptide in the unipept_hapid database has a single species/strain LCA
encoded in the FASTA header as `umgap|<id>|<lca_il>` — DIA-NN surfaces this
as `Protein.Group`, with the trailing `<lca_il>` carrying the NCBI taxon ID.

Spectrum identifier: Run + Precursor.Id (mirrors the hapid module's
build_genome_spectrum_mapping.py).

Output: taxon2spectrum_dic.json  { taxon_id: [spectrum_id, ...] }
"""

import json
import os

import pandas as pd

PARQUET_IN = snakemake.input.parquet
JSON_OUT = snakemake.output[0]
LOG = snakemake.log[0]

os.makedirs(os.path.dirname(JSON_OUT), exist_ok=True)
log = open(LOG, "w")


def logp(msg):
    print(msg)
    print(msg, file=log)


df = pd.read_parquet(PARQUET_IN)
logp(f"Parquet rows: {len(df)}")
if len(df) == 0:
    raise SystemExit(
        f"DIA-NN parquet at {PARQUET_IN} contains 0 rows — "
        "the upstream DIA-NN search produced no peptides. "
        "Inspect the corresponding DIA-NN log under logs/ before re-running."
    )

df = df[df["Proteotypic"] == 1].copy()
logp(f"After Proteotypic == 1 filter: {len(df)} rows")

df["spectrum_id"] = df["Run"].astype(str) + "||" + df["Precursor.Id"].astype(str)
df["taxon_id"] = df["Protein.Group"].astype(str).str.rsplit("|", n=1).str[-1]

taxon2spectrum = (
    df.groupby("taxon_id")["spectrum_id"]
    .apply(lambda s: sorted(set(s)))
    .to_dict()
)
logp(f"Distinct taxa with hits: {len(taxon2spectrum)}")
total_spectra = len(set().union(*taxon2spectrum.values())) if taxon2spectrum else 0
logp(f"Total unique spectra: {total_spectra}")

with open(JSON_OUT, "w") as f:
    json.dump(taxon2spectrum, f)

logp(f"Written to {JSON_OUT}")
log.close()
