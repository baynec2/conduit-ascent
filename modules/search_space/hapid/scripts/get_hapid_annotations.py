#!/usr/bin/env python3
"""
get_hapid_annotations.py

Parse Bakta TSV outputs for selected genomes into mag_annotations.txt format
expected by the mag_annotation module.

Snakemake inputs:
  bakta_dirs — list of Bakta output directories for selected genomes

Snakemake outputs:
  mag_annotations — path to mag_annotations.txt
"""

import os
import re

import pandas as pd

# ── Snakemake bindings ──────────────────────────────────────────────────────
BAKTA_DIRS = list(snakemake.input.bakta_dirs)
OUT_PATH   = snakemake.output.mag_annotations
LOG_PATH   = snakemake.log[0]

os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)

log = open(LOG_PATH, "w")
def logprint(msg):
    print(msg)
    print(msg, file=log, flush=True)

# ── Parse Bakta TSVs ─────────────────────────────────────────────────────────
logprint(f"Parsing Bakta TSVs from {len(BAKTA_DIRS)} directories")


def classify_xref(x):
    if pd.isna(x) or str(x).strip() == "":
        return "unsupported_annotation"
    x = x.strip()
    if x.startswith("SO"):              return "xref_SO"
    if x.startswith("UniRef:UniRef50"): return "xref_uniref50"
    if x.startswith("UniRef:UniRef90"): return "xref_uniref90"
    if x.startswith("UniRef:UniRef100"): return "xref_uniref100"
    if x.startswith("PFAM"):           return "xref_pfam"
    if x.startswith("GO"):             return "xref_go"
    if x.startswith("COG"):            return "xref_cog"
    if x.startswith("EC"):             return "xref_ec"
    return "unsupported_annotation"


def clean_xref_value(x):
    x = re.sub(r"^UniRef:UniRef\d+_", "", x)
    x = re.sub(r"^PFAM:", "", x)
    return x.strip()


frames = []
for d in BAKTA_DIRS:
    mag = os.path.basename(d)
    tsv_path = os.path.join(d, f"{mag}.tsv")
    if not os.path.exists(tsv_path):
        logprint(f"  [{mag}] missing TSV, skipping")
        continue
    try:
        df = pd.read_csv(tsv_path, sep="\t", skiprows=5)
        df["mag"] = mag
        frames.append(df)
        logprint(f"  [{mag}] {len(df)} rows")
    except Exception as e:
        logprint(f"  [{mag}] read error: {e}")

if not frames:
    logprint("No Bakta TSVs found — writing stub")
    stub_cols = ["Locus Tag", "Product", "mag", "xref_uniref100", "xref_uniref90", "xref_uniref50"]
    pd.DataFrame(columns=stub_cols).to_csv(OUT_PATH, sep="\t", index=False)
    log.close()
    raise SystemExit(0)

combined = pd.concat(frames, ignore_index=True)
logprint(f"Total rows: {len(combined)}")

if "DbXrefs" not in combined.columns:
    combined["DbXrefs"] = ""

# Explode DbXrefs — one per row
combined["DbXrefs"] = combined["DbXrefs"].fillna("")
exploded = (
    combined
    .assign(DbXrefs=combined["DbXrefs"].str.split(r",\s*"))
    .explode("DbXrefs")
)
exploded["DbXrefs"] = exploded["DbXrefs"].str.strip()
exploded = exploded[exploded["DbXrefs"] != ""]

# Classify and clean
exploded["xref_name"]  = exploded["DbXrefs"].apply(classify_xref)
exploded["DbXrefs"]    = exploded["DbXrefs"].apply(clean_xref_value)
exploded = exploded[exploded["xref_name"] != "unsupported_annotation"]

# Base table: one row per Locus Tag with metadata columns
id_col = "Locus Tag" if "Locus Tag" in combined.columns else combined.columns[0]
meta_cols = [id_col, "Type", "Name", "Start", "Stop", "Strand",
             "Locus", "Gene", "Product", "mag"]
meta_cols = [c for c in meta_cols if c in combined.columns]

base = combined.drop_duplicates(subset=[id_col])[meta_cols].copy()

# Pivot xref values
if not exploded.empty:
    pivot = (
        exploded
        .groupby([id_col, "xref_name"])["DbXrefs"]
        .apply(lambda vals: ";".join(sorted(set(vals))))
        .unstack("xref_name")
        .reset_index()
    )
    base = base.merge(pivot, on=id_col, how="left")

logprint(f"Output rows: {len(base)}, columns: {list(base.columns)}")
base.to_csv(OUT_PATH, sep="\t", index=False)
logprint(f"Written → {OUT_PATH}")
log.close()
