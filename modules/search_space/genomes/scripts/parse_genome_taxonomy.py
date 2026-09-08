#!/usr/bin/env python3
"""
parse_genome_taxonomy.py

Convert user-provided genome-level taxonomy into conduit's taxonomy.txt format.

Input:  taxonomy.txt — TSV with columns:
          genome  domain  kingdom  phylum  class  order  family  genus  species
Output: genome_taxonomy.txt — one row per genome, adds system-generated columns:
          organism_id (sequential integer, alphabetical order)
          proteome_id, proteome_type, download_info
          genome column is RETAINED for traceability and use by genome_uniprot_headers.py
"""

import os
import pandas as pd

OUT_COLS = [
    "genome", "organism_id", "domain", "kingdom", "phylum", "class", "order",
    "family", "genus", "species", "proteome_id", "proteome_type",
    "download_info"
]

TAXONOMY_COLS = ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"]


def filter_to_selected_genomes(tax_df, selected):
    """Restrict tax_df rows to genomes whose name appears in `selected` (iterable of strings)."""
    return tax_df[tax_df["genome"].astype(str).isin(set(selected))]


def assign_organism_ids(tax_df):
    """Sort by genome name and assign organism_id = 1..N. Returns a new df."""
    out = tax_df.sort_values("genome").reset_index(drop=True)
    out["organism_id"] = range(1, len(out) + 1)
    return out


def fill_taxonomy_defaults(tax_df, logf=None):
    """Add proteome_id/proteome_type/download_info; fill missing taxonomy ranks with 'NA'."""
    out = tax_df.copy()
    out["proteome_id"]   = "NA"
    out["proteome_type"] = "NA"
    out["download_info"] = "user_provided"
    for col in TAXONOMY_COLS:
        if col not in out.columns:
            if logf is not None:
                print(f"WARNING: column '{col}' not found in input — filling with 'NA'", file=logf)
            out[col] = "NA"
    return out


def build_genome_taxonomy(tax_df, selected=None, logf=None):
    """End-to-end pure transform: filter (optional) → assign_organism_ids → fill defaults → OUT_COLS subset."""
    if "genome" not in tax_df.columns:
        raise ValueError("Input taxonomy.txt must have a 'genome' column")
    if selected is not None:
        tax_df = filter_to_selected_genomes(tax_df, selected)
    tax_df = assign_organism_ids(tax_df)
    tax_df = fill_taxonomy_defaults(tax_df, logf=logf)
    return tax_df[OUT_COLS]


def main():
    TAXONOMY_IN   = snakemake.input.taxonomy
    SELECTED_IN   = getattr(snakemake.input, "selected_genomes", None)
    TAXONOMY_OUT  = snakemake.output[0]
    LOG           = snakemake.log[0]

    os.makedirs(os.path.dirname(TAXONOMY_OUT), exist_ok=True)

    log = open(LOG, "w")
    def logprint(msg):
        print(msg)
        print(msg, file=log, flush=True)

    logprint("Reading genome-level taxonomy from taxonomy.txt")
    tax_df = pd.read_csv(TAXONOMY_IN, sep="\t")
    logprint(f"Read {len(tax_df)} genome rows")

    selected = None
    if SELECTED_IN:
        with open(SELECTED_IN) as fh:
            selected = {line.strip() for line in fh if line.strip()}
        before = len(tax_df)
        missing = selected - set(tax_df["genome"].astype(str))
        out_df = build_genome_taxonomy(tax_df, selected=selected, logf=log)
        logprint(
            f"Filtered taxonomy to {len(out_df)}/{before} genomes "
            f"using selected_genomes.txt ({len(selected)} entries)"
        )
        if missing:
            logprint(
                f"WARNING: {len(missing)} selected genomes have no taxonomy row "
                f"(first few: {sorted(missing)[:5]})"
            )
    else:
        out_df = build_genome_taxonomy(tax_df, logf=log)

    out_df.to_csv(TAXONOMY_OUT, sep="\t", index=False)

    logprint(f"Wrote {len(out_df)} genome taxonomy rows → {TAXONOMY_OUT}")
    for _, row in out_df.iterrows():
        logprint(f"  organism_id={row['organism_id']}  genome={row['genome']}  species={row['species']}")

    log.close()


if __name__ == "__main__" or "snakemake" in globals():
    main()
