#!/usr/bin/env python3
"""
parse_hapid_selected_taxonomy.py

Filter genome-level taxonomy to only genomes selected by greedy genome selection,
then convert to conduit's taxonomy.txt format.

Input:
  taxonomy.txt              — TSV with: genome  domain  kingdom  phylum  class  order  family  genus  species
  hapid_greedy_selection.tsv — TSV with: genome  nSpectraCovered  cumulative_pct

Output: hapid_selected_taxonomy.txt — one row per SELECTED genome, with system-generated columns
"""

import os
import pandas as pd

# Snakemake bindings
TAXONOMY_IN  = snakemake.input.taxonomy
SELECTION_IN = snakemake.input.selection
TAXONOMY_OUT = snakemake.output[0]
LOG          = snakemake.log[0]
PERCENT_SPECTRA = snakemake.params.hapid_percent_spectra

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

# ---- Determine selected genomes from greedy selection ----
logprint(f"Reading greedy selection from {SELECTION_IN}")
sel_df = pd.read_csv(SELECTION_IN, sep="\t")
logprint(f"  {len(sel_df)} genomes in greedy ranking")

pct = float(PERCENT_SPECTRA)
above = sel_df[sel_df["cumulative_pct"] >= pct]
cutoff = (above.index[0] + 1) if not above.empty else len(sel_df)
selected_genomes = set(sel_df["genome"].tolist()[:cutoff])
logprint(f"  {len(selected_genomes)} genomes selected at {pct}% spectra threshold")

# ---- Read and filter taxonomy ----
logprint(f"Reading genome-level taxonomy from {TAXONOMY_IN}")
tax_df = pd.read_csv(TAXONOMY_IN, sep="\t")
logprint(f"  {len(tax_df)} total genome rows")

if "genome" not in tax_df.columns:
    raise ValueError("Input taxonomy.txt must have a 'genome' column")

tax_df = tax_df[tax_df["genome"].isin(selected_genomes)]
logprint(f"  {len(tax_df)} genomes after filtering to selected")

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

logprint(f"Wrote {len(out_df)} selected genome taxonomy rows → {TAXONOMY_OUT}")
for _, row in out_df.iterrows():
    logprint(f"  organism_id={row['organism_id']}  genome={row['genome']}  species={row['species']}")

log.close()
