#!/usr/bin/env python3

from Bio import SeqIO
import pandas as pd
import os

# -----------------------------
# Snakemake bindings
# -----------------------------
BAKTA_DIRS    = snakemake.input.bakta_dirs
TAXONOMY_PATH = snakemake.params.taxonomy   # mag_taxonomy.txt (has genome + organism_id cols)

FASTA_OUT = snakemake.output.fasta

LOG = snakemake.log[0]

os.makedirs(os.path.dirname(FASTA_OUT), exist_ok=True)

log = open(LOG, "w")
def logprint(msg):
    print(msg)
    print(msg, file=log)

# -----------------------------
# Map MAG name -> its bakta output dir, directly from the snakemake input.
# Avoids reconstructing paths from a derived BASE_DIR, which broke for the
# 1-MAG case (commonpath returns the dir itself, not its parent) and for the
# 0-MAG case (IndexError on BAKTA_DIRS[0]).
# -----------------------------
MAG_DIR_BY_NAME = {os.path.basename(d): d for d in BAKTA_DIRS}

# -----------------------------
# Load taxonomy lookup (genome → species_name, organism_id)
# -----------------------------
tax_df = pd.read_csv(TAXONOMY_PATH, sep="\t")
GENOME_SPECIES   = dict(zip(tax_df["genome"].astype(str), tax_df["species"].astype(str)))
GENOME_ORG_ID    = dict(zip(tax_df["genome"].astype(str), tax_df["organism_id"].astype(str)))
logprint(f"Loaded taxonomy for {len(GENOME_SPECIES)} genomes from {TAXONOMY_PATH}")

# -----------------------------
# Helper functions
# -----------------------------

def preprocess_bakta_annotations(MAG_dir, MAG):
    """
    Extract CDS-only fastas and metadata from bakta output
    """
    faa_path = os.path.join(MAG_dir, f"{MAG}.faa")
    faa_clean_path = os.path.join(MAG_dir, f"{MAG}_cds.faa")

    tsv_path = os.path.join(MAG_dir, f"{MAG}.tsv")
    cds_tsv_path = os.path.join(MAG_dir, f"{MAG}_cds.tsv")

    if not os.path.exists(tsv_path):
        logprint(f"[{MAG}] Missing TSV: {tsv_path}")
        return None, None

    if not os.path.exists(faa_path):
        logprint(f"[{MAG}] Missing FAA: {faa_path}")
        return None, None

    # Bakta TSV has 5-line header
    meta = pd.read_csv(tsv_path, sep="\t", skiprows=5)
    cds_meta = meta[meta["Type"] == "cds"]

    records = list(SeqIO.parse(faa_path, "fasta"))

    faa_ids = {r.id for r in records}
    cds_ids = set(cds_meta["Locus Tag"])

    keep_ids = faa_ids.intersection(cds_ids)
    logprint(f"[{MAG}] Found {len(keep_ids)} CDS IDs")

    filtered_records = [r for r in records if r.id in keep_ids]
    cds_meta_filtered = cds_meta[cds_meta["Locus Tag"].isin(keep_ids)]

    SeqIO.write(filtered_records, faa_clean_path, "fasta")
    cds_meta_filtered.to_csv(cds_tsv_path, sep="\t", index=False)

    return faa_clean_path, cds_tsv_path


def format_uniprot_metadata(MAG_dir, MAG):
    """
    Creates a dataframe with UniProt-style metadata. organism_id is looked up
    from the taxonomy file (not computed from the genome name).
    """
    cds_tsv_path = os.path.join(MAG_dir, f"{MAG}_cds.tsv")

    if not os.path.exists(cds_tsv_path):
        logprint(f"[{MAG}] Missing CDS TSV: {cds_tsv_path}")
        return None

    cds_meta = pd.read_csv(cds_tsv_path, sep="\t")
    locus_tags = cds_meta["Locus Tag"].tolist()

    df = pd.DataFrame(index=locus_tags)
    df["Product"] = cds_meta.set_index("Locus Tag")["Product"].reindex(locus_tags).fillna("NA")
    df["Gene"]    = cds_meta.set_index("Locus Tag")["Gene"].reindex(locus_tags).fillna("NA")

    species_name = GENOME_SPECIES.get(MAG, "Unknown species")
    org_id       = GENOME_ORG_ID.get(MAG, "0")

    # MGnify metadata has empty species cells for some MAGs -> pandas reads NaN.
    if pd.isna(species_name) or not isinstance(species_name, str):
        species_name = "Unknown species"

    parts = species_name.split(" ")
    genus   = parts[0] if len(parts) > 0 else "UNK"
    species = parts[1] if len(parts) > 1 else "UNK"

    genus_abbr   = genus[:3].upper() if genus else "UNK"
    species_abbr = species[:2].upper() if species else "UNK"

    df["Species_Locus_Tag"]  = [f"{lt}_{genus_abbr}{species_abbr}" for lt in locus_tags]
    df["Species_Name"]       = species_name
    df["Organism_Identifier"] = org_id  # from taxonomy lookup

    return df


def replace_faa_headers(MAG_dir, MAG, df):
    """
    Replace FASTA headers with UniProt-style annotation.
    """
    if df is None:
        return None

    faa_in  = os.path.join(MAG_dir, f"{MAG}_cds.faa")
    faa_out = os.path.join(MAG_dir, f"{MAG}_uniprot.faa")

    if not os.path.exists(faa_in):
        logprint(f"[{MAG}] Missing CDS FAA: {faa_in}")
        return None

    records = list(SeqIO.parse(faa_in, "fasta"))
    updated = 0

    for r in records:
        if r.id in df.index:
            meta = df.loc[r.id]
            header = (
                f"tr|{r.id}|{meta['Species_Locus_Tag']} "
                f"{meta['Product']} OS={meta['Species_Name']} "
                f"OX={meta['Organism_Identifier']} GN={meta['Gene']}"
            )
            r.id = header
            r.description = ""
            updated += 1

    SeqIO.write(records, faa_out, "fasta")
    logprint(f"[{MAG}] Updated {updated} headers → {faa_out}")

    return faa_out


def concatenate_uniprot_fastas(mag_dirs_by_name, output_path):
    """
    Combines all MAG uniprot fastas into one file. mag_dirs_by_name is the
    {basename: full_path} mapping of bakta output dirs.
    """
    files = []
    for MAG, d in mag_dirs_by_name.items():
        f = os.path.join(d, f"{MAG}_uniprot.faa")
        if os.path.exists(f):
            files.append(f)

    if not files:
        logprint("No uniprot FASTAs found — writing empty combined FASTA.")
        open(output_path, "w").close()
        return

    with open(output_path, "w") as out:
        for f in files:
            with open(f) as inp:
                out.write(inp.read())

    logprint(f"Wrote combined FASTA → {output_path}")


# -----------------------------
# MAIN
# -----------------------------
def main():
    MAG_list = list(GENOME_SPECIES.keys())

    for MAG in MAG_list:
        MAG_dir = MAG_DIR_BY_NAME.get(MAG)
        if not MAG_dir or not os.path.isdir(MAG_dir):
            logprint(f"[{MAG}] Skipping (no bakta directory in inputs)")
            continue

        logprint(f"Processing {MAG}")

        preprocess_bakta_annotations(MAG_dir, MAG)
        df = format_uniprot_metadata(MAG_dir, MAG)
        replace_faa_headers(MAG_dir, MAG, df)

    concatenate_uniprot_fastas(MAG_DIR_BY_NAME, FASTA_OUT)

    logprint("MAG UniProt processing complete.")
    log.close()


main()
