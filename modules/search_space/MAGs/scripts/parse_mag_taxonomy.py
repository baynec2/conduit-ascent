#!/usr/bin/env python3
"""
parse_mag_taxonomy.py

Convert user-provided genome-level taxonomy into conduit's taxonomy.txt format.

Input:  taxonomy.txt — TSV with columns:
          genome  domain  kingdom  phylum  class  order  family  genus  species
Output: mag_taxonomy.txt — one row per genome, adds system-generated columns:
          organism_id (sequential integer, alphabetical order)
          proteome_id, proteome_type, download_info, organism_type
          genome column is RETAINED for traceability and use by MAG_uniprot_headers.py
"""

import os
import pandas as pd

# Snakemake bindings
TAXONOMY_IN  = snakemake.input.taxonomy
TAXONOMY_OUT = snakemake.output[0]
LOG          = snakemake.log[0]

os.makedirs(os.path.dirname(TAXONOMY_OUT), exist_ok=True)

log = open(LOG, "w")
def logprint(msg):
    print(msg)
    print(msg, file=log, flush=True)

# Columns for output (genome col retained for traceability)
OUT_COLS = [
    "genome", "organism_id", "domain", "kingdom", "phylum", "class", "order",
    "family", "genus", "species", "proteome_id", "proteome_type",
    "download_info", "organism_type"
]

TAXONOMY_COLS = ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"]

logprint("Reading genome-level taxonomy from taxonomy.txt")

tax_df = pd.read_csv(TAXONOMY_IN, sep="\t")
logprint(f"Read {len(tax_df)} genome rows")

if "genome" not in tax_df.columns:
    raise ValueError("Input taxonomy.txt must have a 'genome' column")

# Assign organism_id as sequential integers by alphabetical genome order
tax_df = tax_df.sort_values("genome").reset_index(drop=True)
tax_df["organism_id"] = range(1, len(tax_df) + 1)

# Fill in system-generated columns
tax_df["proteome_id"]   = "NA"
tax_df["proteome_type"] = "NA"
tax_df["download_info"] = "user_provided"
tax_df["organism_type"] = "microbiome"

# Ensure all expected taxonomy columns exist (fill NA if absent)
for col in TAXONOMY_COLS:
    if col not in tax_df.columns:
        logprint(f"WARNING: column '{col}' not found in input — filling with 'NA'")
        tax_df[col] = "NA"

out_df = tax_df[OUT_COLS]
out_df.to_csv(TAXONOMY_OUT, sep="\t", index=False)

logprint(f"Wrote {len(out_df)} genome taxonomy rows → {TAXONOMY_OUT}")
for _, row in out_df.iterrows():
    logprint(f"  organism_id={row['organism_id']}  genome={row['genome']}  species={row['species']}")

log.close()
