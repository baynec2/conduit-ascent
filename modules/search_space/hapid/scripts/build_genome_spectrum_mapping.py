#!/usr/bin/env python3
"""
build_genome_spectrum_mapping.py

Map DIA-NN profiling search results (parquet) to genomes using:
  1. protein2genome_dic.json  ({protein_id: genome_id})
  2. CD-HIT cluster file (.clstr) to expand cluster members (skipped in
     precompiled mode when no .clstr exists)

DIA-NN parquet columns used:
  - Run + Precursor.Id  → unique spectrum identifier
  - Protein.Ids         → semicolon-separated protein IDs

Output: genome2spectrum_dic.json  { genome_id: [spectrum_id, ...] }
"""

import json
import os
import sys
import pandas as pd

# Snakemake bindings
PARQUET_IN         = snakemake.input.parquet
PROTEIN_GENOME_DIC = snakemake.input.protein2genome
CLSTR_IN           = snakemake.input.get("clstr", None)
JSON_OUT           = snakemake.output[0]
LOG                = snakemake.log[0]

os.makedirs(os.path.dirname(JSON_OUT), exist_ok=True)

log = open(LOG, "w")
def logprint(msg):
    print(msg)
    print(msg, file=log)


def parse_cdhit_clstr(clstr_path):
    """
    Parse a CD-HIT .clstr file. Returns dict:
      { representative_protein_id: [all_member_protein_ids] }
    All IDs are extracted by stripping '>....' and '...' prefixes.
    """
    repr2members = {}
    current_members = []
    current_repr = None

    with open(clstr_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                if current_repr is not None:
                    repr2members[current_repr] = current_members
                current_members = []
                current_repr = None
            else:
                # e.g.: 0    123aa, >genome|prot... *
                parts = line.split()
                if len(parts) < 3:
                    continue
                # extract protein id from ">genome|prot..."
                raw_id = parts[2].lstrip(">").rstrip(".")
                # strip genome prefix (genome_id|protein_id → protein_id)
                if "|" in raw_id:
                    prot_id = raw_id.split("|", 1)[1]
                else:
                    prot_id = raw_id
                current_members.append(prot_id)
                if line.endswith("*"):
                    current_repr = prot_id
        if current_repr is not None:
            repr2members[current_repr] = current_members

    return repr2members


# Load protein→genome dictionary
with open(PROTEIN_GENOME_DIC) as f:
    protein2genome = json.load(f)

logprint(f"Loaded {len(protein2genome)} protein→genome mappings")

# Build CD-HIT cluster expansion map (repr → all members)
repr2members = {}
if CLSTR_IN and os.path.exists(CLSTR_IN):
    repr2members = parse_cdhit_clstr(CLSTR_IN)
    logprint(f"Loaded {len(repr2members)} CD-HIT clusters for expansion")
else:
    logprint("No CD-HIT cluster file — skipping cluster expansion")

# Read DIA-NN parquet
df = pd.read_parquet(PARQUET_IN)
logprint(f"Parquet rows: {len(df)}")

# Build unique spectrum ID from Run + Precursor.Id
df["spectrum_id"] = df["Run"].astype(str) + "||" + df["Precursor.Id"].astype(str)

# Build protein2spectrum mapping
protein2spectrum = {}  # { protein_id: [spectrum_id, ...] }
for _, row in df.iterrows():
    spec_id = row["spectrum_id"]
    protein_ids_str = row["Protein.Ids"]
    if not isinstance(protein_ids_str, str):
        continue
    for raw_prot in protein_ids_str.split(";"):
        raw_prot = raw_prot.strip()
        # DIA-NN headers are genome_id|protein_id — extract protein_id
        if "|" in raw_prot:
            prot_id = raw_prot.split("|", 1)[1]
        else:
            prot_id = raw_prot
        if prot_id not in protein2spectrum:
            protein2spectrum[prot_id] = []
        protein2spectrum[prot_id].append(spec_id)

# Deduplicate
protein2spectrum = {k: list(set(v)) for k, v in protein2spectrum.items()}
logprint(f"Unique proteins with spectrum hits: {len(protein2spectrum)}")

# Expand via CD-HIT clusters: if a representative has hits, credit all members
protein2spectrum_expanded = dict(protein2spectrum)
for repr_prot, members in repr2members.items():
    if repr_prot in protein2spectrum:
        repr_spectra = protein2spectrum[repr_prot]
        for member in members:
            if member != repr_prot:
                if member not in protein2spectrum_expanded:
                    protein2spectrum_expanded[member] = []
                protein2spectrum_expanded[member] = list(
                    set(protein2spectrum_expanded[member]) | set(repr_spectra)
                )

logprint(f"After cluster expansion: {len(protein2spectrum_expanded)} proteins with spectra")

# Aggregate to genome level
genome2spectrum = {}
for prot_id, spectra in protein2spectrum_expanded.items():
    genome_id = protein2genome.get(prot_id)
    if genome_id is None:
        continue
    if genome_id not in genome2spectrum:
        genome2spectrum[genome_id] = []
    genome2spectrum[genome_id].extend(spectra)

genome2spectrum = {k: list(set(v)) for k, v in genome2spectrum.items()}
logprint(f"Genome→spectrum mapping: {len(genome2spectrum)} genomes")
for gid, spectra in sorted(genome2spectrum.items(), key=lambda x: -len(x[1]))[:10]:
    logprint(f"  {gid}: {len(spectra)} spectra")

with open(JSON_OUT, "w") as f:
    json.dump(genome2spectrum, f)

logprint(f"Written to {JSON_OUT}")
log.close()
