#!/usr/bin/env python3
"""
map_families_to_species.py

Given detected family names from the first-pass FDR inference, look up the
user-provided taxonomy to find all species belonging to those families.

Replaces the taxonkit-based rule in the original peptidotyping module.

Input:
  - detected_family_taxa_ids.txt (TSV: ncbi_taxonomy_id, detected_taxonomy)
  - taxonomy.txt (user-provided genome taxonomy)

Output:
  - families_to_species.txt (one underscore-joined species name per line)
"""

import os

import pandas as pd

# Snakemake bindings
DETECTED_FAMILIES_PATH = snakemake.input.detected_families
TAXONOMY_PATH = snakemake.input.taxonomy
OUTPUT_PATH = snakemake.output.families_to_species
LOG_PATH = snakemake.log[0]

os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)

log = open(LOG_PATH, "w")


def logprint(msg):
    print(msg)
    print(msg, file=log, flush=True)


def underscore_name(name):
    if pd.isna(name) or str(name).strip().upper() == "NA":
        return "NA"
    return str(name).strip().replace(" ", "_")


def main():
    # Read detected families (column 1 = family identifier, underscore-joined)
    det_df = pd.read_csv(DETECTED_FAMILIES_PATH, sep="\t")
    if len(det_df) == 0:
        logprint("No families detected — writing empty species list")
        open(OUTPUT_PATH, "w").close()
        log.close()
        return

    detected_families = set(det_df.iloc[:, 0].astype(str).tolist())
    logprint(f"Detected {len(detected_families)} families: {detected_families}")

    # Read user taxonomy
    tax_df = pd.read_csv(TAXONOMY_PATH, sep="\t")
    logprint(f"Taxonomy has {len(tax_df)} genomes")

    # Find all species in detected families
    species_set = set()
    for _, row in tax_df.iterrows():
        family_u = underscore_name(row.get("family", "NA"))
        species_u = underscore_name(row.get("species", "NA"))

        if family_u in detected_families and species_u != "NA":
            species_set.add(species_u)

    logprint(f"Found {len(species_set)} species in detected families")

    with open(OUTPUT_PATH, "w") as f:
        for sp in sorted(species_set):
            f.write(f"{sp}\n")

    logprint("map_families_to_species.py complete")
    log.close()


main()
