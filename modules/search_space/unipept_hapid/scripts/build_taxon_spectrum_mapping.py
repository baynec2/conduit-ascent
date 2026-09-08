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


def build_taxon_spectrum_mapping(df):
    """Pure transform: filter proteotypic precursors, group spectrum_ids by taxon.

    Expects columns: Proteotypic, Run, Precursor.Id, Protein.Group.
    Returns dict {taxon_id: sorted list of unique spectrum_ids}.
    """
    df = df[df["Proteotypic"] == 1].copy()
    df["spectrum_id"] = df["Run"].astype(str) + "||" + df["Precursor.Id"].astype(str)
    df["taxon_id"] = df["Protein.Group"].astype(str).str.rsplit("|", n=1).str[-1]
    return (
        df.groupby("taxon_id")["spectrum_id"]
          .apply(lambda s: sorted(set(s)))
          .to_dict()
    )


def main():
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
        # Genuine no-detection: the HAPiID first pass identified no peptides for
        # this sample. Rather than hard-fail, emit an empty taxon map so the run
        # resolves to an empty (no-detection) conduit downstream — see the
        # database.fasta emptiness checkpoint in database_processing. Inspect the
        # DIA-NN first-pass log under logs/ if a non-empty result was expected.
        logp(
            f"DIA-NN parquet at {PARQUET_IN} contains 0 rows — no peptides "
            "identified. Writing an empty taxon map (run resolves to an empty "
            "conduit)."
        )
        with open(JSON_OUT, "w") as f:
            json.dump({}, f)
        logp(f"Written empty map to {JSON_OUT}")
        log.close()
        return

    taxon2spectrum = build_taxon_spectrum_mapping(df)
    logp(f"After Proteotypic == 1 filter: {sum(len(v) for v in taxon2spectrum.values())} precursor-rows kept across {len(taxon2spectrum)} taxa")
    total_spectra = len(set().union(*taxon2spectrum.values())) if taxon2spectrum else 0
    logp(f"Total unique spectra: {total_spectra}")

    with open(JSON_OUT, "w") as f:
        json.dump(taxon2spectrum, f)

    logp(f"Written to {JSON_OUT}")
    log.close()


if __name__ == "__main__" or "snakemake" in globals():
    main()
