#!/usr/bin/env python3
"""
parse_mgnify_metadata.py

Parse MGnify genomes-all_metadata.tsv to:
  1. Filter to species-representative genomes (Genome == Species_rep)
  2. Optionally filter by GTDB lineage substring
  3. Optionally cap number of genomes
  4. Convert GTDB taxonomy to conduit's taxonomy.txt format
  5. Write species_representatives.txt (one accession per line)

GTDB lineage format:
  d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Blautia_A;s__Blautia_A faecis
"""

import os
import pandas as pd

# Snakemake bindings
METADATA_IN       = snakemake.input.metadata
TAXONOMY_OUT      = snakemake.output.taxonomy
REPRESENTATIVES   = snakemake.output.representatives
LOG               = snakemake.log[0]
TAXONOMY_FILTER   = snakemake.params.taxonomy_filter
MAX_GENOMES       = int(snakemake.params.max_genomes)

os.makedirs(os.path.dirname(TAXONOMY_OUT), exist_ok=True)
os.makedirs(os.path.dirname(REPRESENTATIVES), exist_ok=True)

log = open(LOG, "w")
def logprint(msg):
    print(msg)
    print(msg, file=log, flush=True)

# GTDB rank prefixes in order
GTDB_PREFIXES = ["d__", "p__", "c__", "o__", "f__", "g__", "s__"]
TAXONOMY_COLS = ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"]


def parse_gtdb_lineage(lineage):
    """Parse a GTDB lineage string into a dict of taxonomy columns."""
    result = {col: "NA" for col in TAXONOMY_COLS}

    if pd.isna(lineage) or not lineage.strip():
        return result

    parts = lineage.split(";")
    rank_map = {}
    for part in parts:
        part = part.strip()
        for prefix in GTDB_PREFIXES:
            if part.startswith(prefix):
                rank_map[prefix] = part[len(prefix):].strip() or "NA"
                break

    result["domain"]  = rank_map.get("d__", "NA")
    result["kingdom"] = result["domain"]  # No kingdom in GTDB; mirror domain
    result["phylum"]  = rank_map.get("p__", "NA")
    result["class"]   = rank_map.get("c__", "NA")
    result["order"]   = rank_map.get("o__", "NA")
    result["family"]  = rank_map.get("f__", "NA")
    result["genus"]   = rank_map.get("g__", "NA")
    result["species"] = rank_map.get("s__", "NA")

    return result


logprint("Reading MGnify metadata")
df = pd.read_csv(METADATA_IN, sep="\t", low_memory=False)
logprint(f"Total genomes in metadata: {len(df)}")

# Validate required columns
for col in ["Genome", "Species_rep", "Lineage"]:
    if col not in df.columns:
        raise ValueError(
            f"Expected column '{col}' not found in metadata. "
            f"Available columns: {list(df.columns)}"
        )

# Filter to species representatives
species_reps = df[df["Genome"] == df["Species_rep"]].copy()
logprint(f"Species representatives: {len(species_reps)}")

# Apply optional taxonomy filter
if TAXONOMY_FILTER and TAXONOMY_FILTER != "FALSE" and TAXONOMY_FILTER is not False:
    before = len(species_reps)
    species_reps = species_reps[
        species_reps["Lineage"].str.contains(TAXONOMY_FILTER, na=False)
    ]
    logprint(f"After taxonomy filter '{TAXONOMY_FILTER}': {len(species_reps)} (removed {before - len(species_reps)})")

# Apply optional max genomes cap
if MAX_GENOMES > 0 and len(species_reps) > MAX_GENOMES:
    species_reps = species_reps.head(MAX_GENOMES)
    logprint(f"Capped to {MAX_GENOMES} genomes")

if len(species_reps) == 0:
    raise ValueError("No species representatives remaining after filtering. Check mgnify_taxonomy_filter.")

# Parse GTDB lineage into taxonomy columns
logprint("Parsing GTDB lineage into taxonomy columns")
tax_records = []
for _, row in species_reps.iterrows():
    parsed = parse_gtdb_lineage(row["Lineage"])
    parsed["genome"] = row["Genome"]
    tax_records.append(parsed)

tax_df = pd.DataFrame(tax_records)
tax_df = tax_df.sort_values("genome").reset_index(drop=True)

# Write taxonomy.txt
tax_df = tax_df[["genome"] + TAXONOMY_COLS]
tax_df.to_csv(TAXONOMY_OUT, sep="\t", index=False)
logprint(f"Wrote taxonomy.txt with {len(tax_df)} genomes -> {TAXONOMY_OUT}")

# Write species_representatives.txt
accessions = tax_df["genome"].tolist()
with open(REPRESENTATIVES, "w") as f:
    for acc in accessions:
        f.write(acc + "\n")
logprint(f"Wrote species_representatives.txt with {len(accessions)} accessions -> {REPRESENTATIVES}")

# Log summary
for _, row in tax_df.head(10).iterrows():
    logprint(f"  {row['genome']}: {row['species']} ({row['phylum']})")
if len(tax_df) > 10:
    logprint(f"  ... and {len(tax_df) - 10} more")

log.close()
