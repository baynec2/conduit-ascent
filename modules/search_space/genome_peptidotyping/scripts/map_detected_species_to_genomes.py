#!/usr/bin/env python3
"""
map_detected_species_to_genomes.py

Maps detected species/strain names from the second-pass FDR inference back to
input genome names using the user-provided taxonomy.

This is the checkpoint script — its output controls which genomes the
downstream MAGs rules process.

Input:
  - detected_species_strain_taxa_ids.txt (TSV: ncbi_taxonomy_id, detected_taxonomy)
  - taxonomy.txt (user-provided genome taxonomy)

Output:
  - detected_genomes.txt (one genome name per line)
"""

import os

import pandas as pd

# Snakemake bindings
DETECTED_SPECIES_PATH = snakemake.input.detected_species
TAXONOMY_PATH = snakemake.input.taxonomy
OUTPUT_PATH = snakemake.output.detected_genomes
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
    # Read detected species (column 1 = species identifier, underscore-joined)
    det_df = pd.read_csv(DETECTED_SPECIES_PATH, sep="\t")
    if len(det_df) == 0:
        logprint("WARNING: No species detected — writing empty genome list")
        open(OUTPUT_PATH, "w").close()
        log.close()
        return

    detected_species = set(det_df.iloc[:, 0].astype(str).tolist())
    logprint(f"Detected {len(detected_species)} species/strains")

    # Read user taxonomy
    tax_df = pd.read_csv(TAXONOMY_PATH, sep="\t")
    logprint(f"Taxonomy has {len(tax_df)} genomes")

    # Match detected species to genomes
    selected_genomes = []
    for _, row in tax_df.iterrows():
        genome = str(row["genome"])
        species_u = underscore_name(row.get("species", "NA"))

        if species_u in detected_species:
            selected_genomes.append(genome)

    logprint(f"Selected {len(selected_genomes)} genomes from {len(tax_df)} total")

    with open(OUTPUT_PATH, "w") as f:
        for g in sorted(selected_genomes):
            f.write(f"{g}\n")

    if len(selected_genomes) == 0:
        logprint("WARNING: No genomes matched any detected species")
    else:
        for g in sorted(selected_genomes):
            logprint(f"  Selected: {g}")

    logprint("map_detected_species_to_genomes.py complete")
    log.close()


main()
