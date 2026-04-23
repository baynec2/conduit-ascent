#!/usr/bin/env python3
"""
hapid_uniprot_headers.py

Convert Bakta annotation output for selected HAPiID genomes into UniProt-style
FASTA headers. Adapted from MAGs/scripts/MAG_uniprot_headers.py.

Key differences from MAG_uniprot_headers.py:
  - Species name and organism_id are both read from hapid_taxonomy.txt
    (the processed taxonomy output of parse_hapid_taxonomy.py).
  - organism_id is a sequential integer assigned by parse_hapid_taxonomy.py —
    never computed independently here.
  - "mag" terminology replaced by "genome"
"""

from Bio import SeqIO
import pandas as pd
import os
import re
import glob

# Snakemake bindings
BAKTA_DIRS     = snakemake.input.bakta_dirs
HAPID_TAXONOMY = snakemake.params.hapid_taxonomy

FASTA_OUT = snakemake.output.fasta
GO_OUT    = snakemake.output.go
KEGG_OUT  = snakemake.output.kegg

LOG = snakemake.log[0]

os.makedirs(os.path.dirname(FASTA_OUT), exist_ok=True)

log = open(LOG, "w")
def logprint(msg):
    print(msg)
    print(msg, file=log)


# Build genome_id → species_name and genome_id → organism_id from hapid_taxonomy.txt
def build_genome_maps(hapid_taxonomy_path):
    tax_df = pd.read_csv(hapid_taxonomy_path, sep="\t")
    genome_species   = dict(zip(tax_df["genome"].astype(str), tax_df["species"].astype(str)))
    genome_org_id    = dict(zip(tax_df["genome"].astype(str), tax_df["organism_id"].astype(str)))
    return genome_species, genome_org_id


GENOME_SPECIES, GENOME_ORGANISM_ID = build_genome_maps(HAPID_TAXONOMY)
logprint(f"Loaded taxonomy for {len(GENOME_SPECIES)} genomes from {HAPID_TAXONOMY}")

BASE_DIR = os.path.dirname(BAKTA_DIRS[0])


def preprocess_bakta_annotations(genome_dir, genome):
    faa_path      = os.path.join(genome_dir, f"{genome}.faa")
    faa_clean     = os.path.join(genome_dir, f"{genome}_cds.faa")
    tsv_path      = os.path.join(genome_dir, f"{genome}.tsv")
    cds_tsv_path  = os.path.join(genome_dir, f"{genome}_cds.tsv")

    if not os.path.exists(tsv_path):
        logprint(f"[{genome}] Missing TSV: {tsv_path}")
        return None, None
    if not os.path.exists(faa_path):
        logprint(f"[{genome}] Missing FAA: {faa_path}")
        return None, None

    meta     = pd.read_csv(tsv_path, sep="\t", skiprows=5)
    cds_meta = meta[meta["Type"] == "cds"]
    records  = list(SeqIO.parse(faa_path, "fasta"))

    faa_ids  = {r.id for r in records}
    cds_ids  = set(cds_meta["Locus Tag"])
    keep_ids = faa_ids.intersection(cds_ids)
    logprint(f"[{genome}] {len(keep_ids)} CDS IDs kept")

    filtered = [r for r in records if r.id in keep_ids]
    cds_meta[cds_meta["Locus Tag"].isin(keep_ids)].to_csv(cds_tsv_path, sep="\t", index=False)
    SeqIO.write(filtered, faa_clean, "fasta")

    return faa_clean, cds_tsv_path


def extract_go_kegg(dbxref_string):
    if pd.isna(dbxref_string):
        return [], []
    gos   = re.findall(r"GO:\d+", dbxref_string)
    keggs = re.findall(r"(?:KEGG:)?(K\d+)", dbxref_string)
    return gos, keggs


def format_uniprot_metadata(genome_dir, genome):
    cds_tsv = os.path.join(genome_dir, f"{genome}_cds.tsv")
    if not os.path.exists(cds_tsv):
        logprint(f"[{genome}] Missing CDS TSV")
        return None

    cds_meta   = pd.read_csv(cds_tsv, sep="\t")
    locus_tags = cds_meta["Locus Tag"].tolist()

    df = pd.DataFrame(index=locus_tags)
    df["Product"] = cds_meta.set_index("Locus Tag")["Product"].reindex(locus_tags).fillna("NA")
    df["Gene"]    = cds_meta.set_index("Locus Tag")["Gene"].reindex(locus_tags).fillna("NA")

    species_name = GENOME_SPECIES.get(genome, "Unknown species")
    org_id       = GENOME_ORGANISM_ID.get(genome, "0")

    parts        = species_name.split(" ")
    genus        = parts[0] if len(parts) > 0 else "UNK"
    species_part = parts[1] if len(parts) > 1 else "UNK"

    genus_abbr   = genus[:3].upper()
    species_abbr = species_part[:2].upper()

    df["Species_Locus_Tag"]   = [f"{lt}_{genus_abbr}{species_abbr}" for lt in locus_tags]
    df["Species_Name"]        = species_name
    df["Organism_Identifier"] = org_id  # from taxonomy lookup, not computed

    go_terms, kegg_terms = [], []
    if "DbXrefs" in cds_meta.columns:
        for x in cds_meta["DbXrefs"]:
            g, k = extract_go_kegg(x)
            go_terms.append(g)
            kegg_terms.append(k)
    else:
        go_terms   = [[] for _ in locus_tags]
        kegg_terms = [[] for _ in locus_tags]

    df["GO_Terms"]   = go_terms
    df["KEGG_Terms"] = kegg_terms

    return df


def replace_faa_headers(genome_dir, genome, df):
    if df is None:
        return None

    faa_in  = os.path.join(genome_dir, f"{genome}_cds.faa")
    faa_out = os.path.join(genome_dir, f"{genome}_uniprot.faa")

    if not os.path.exists(faa_in):
        logprint(f"[{genome}] Missing CDS FAA")
        return None

    records = list(SeqIO.parse(faa_in, "fasta"))
    updated = 0
    for r in records:
        if r.id in df.index:
            meta   = df.loc[r.id]
            header = (
                f"tr|{r.id}|{meta['Species_Locus_Tag']} "
                f"{meta['Product']} OS={meta['Species_Name']} "
                f"OX={meta['Organism_Identifier']} GN={meta['Gene']}"
            )
            r.id          = header
            r.description = ""
            updated += 1

    SeqIO.write(records, faa_out, "fasta")
    logprint(f"[{genome}] Updated {updated} headers → {faa_out}")
    return faa_out


def concatenate_uniprot_fastas(base_dir, output_path):
    files = []
    for genome in os.listdir(base_dir):
        f = os.path.join(base_dir, genome, f"{genome}_uniprot.faa")
        if os.path.exists(f):
            files.append(f)

    if not files:
        logprint("No uniprot FASTAs found.")
        return

    with open(output_path, "w") as out:
        for f in files:
            with open(f) as inp:
                out.write(inp.read())

    logprint(f"Wrote combined FASTA → {output_path}")


def write_go_kegg_annotation_files(df):
    df[["GO_Terms"]].to_csv(GO_OUT, sep="\t", header=False)
    df[["KEGG_Terms"]].to_csv(KEGG_OUT, sep="\t", header=False)
    logprint(f"Wrote GO annotations → {GO_OUT}")
    logprint(f"Wrote KEGG annotations → {KEGG_OUT}")


# MAIN
all_meta = []
for bakta_dir in BAKTA_DIRS:
    genome = os.path.basename(bakta_dir)
    logprint(f"Processing {genome}")

    preprocess_bakta_annotations(bakta_dir, genome)
    df = format_uniprot_metadata(bakta_dir, genome)
    replace_faa_headers(bakta_dir, genome, df)

    if df is not None:
        df["genome"] = genome
        all_meta.append(df)

if not all_meta:
    logprint("No genome metadata found — nothing to export.")
    log.close()
    raise RuntimeError("No genomes were processed")

full_df = pd.concat(all_meta)
write_go_kegg_annotation_files(full_df)
concatenate_uniprot_fastas(BASE_DIR, FASTA_OUT)

logprint("HAPiID UniProt + GO/KEGG processing complete.")
log.close()
