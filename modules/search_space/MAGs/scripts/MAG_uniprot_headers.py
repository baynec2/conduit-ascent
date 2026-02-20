#!/usr/bin/env python3

from Bio import SeqIO
import pandas as pd
import os
import re

# -----------------------------
# Snakemake bindings
# -----------------------------
BAKTA_DIRS = snakemake.input.bakta_dirs
MAG_METADATA_PATH = snakemake.params.mag_metadata

FASTA_OUT = snakemake.output.fasta
GO_OUT = snakemake.output.go
KEGG_OUT = snakemake.output.kegg

LOG = snakemake.log[0]

# Create directories for outputs
os.makedirs(os.path.dirname(FASTA_OUT), exist_ok=True)
os.makedirs(os.path.dirname(GO_OUT), exist_ok=True)
os.makedirs(os.path.dirname(KEGG_OUT), exist_ok=True)

log = open(LOG, "w")
def logprint(msg):
    print(msg)
    print(msg, file=log)

# -----------------------------
# Determine base bakta directory
# -----------------------------
BASE_DIR = os.path.commonpath(BAKTA_DIRS)


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



def extract_go_kegg(dbxref_string):
    """
    Extract GO:xxxxxxx and KEGG K numbers.
    """
    if pd.isna(dbxref_string):
        return [], []

    gos = re.findall(r"GO:\d+", dbxref_string)
    keggs = re.findall(r"(?:KEGG:)?(K\d+)", dbxref_string)

    return gos, keggs



def format_uniprot_metadata(MAG_dir, MAG, MAG_info):
    """
    Creates a dataframe with UniProt-style metadata and GO/KEGG lists.
    """
    cds_tsv_path = os.path.join(MAG_dir, f"{MAG}_cds.tsv")

    if not os.path.exists(cds_tsv_path):
        logprint(f"[{MAG}] Missing CDS TSV: {cds_tsv_path}")
        return None

    cds_meta = pd.read_csv(cds_tsv_path, sep="\t")
    locus_tags = cds_meta["Locus Tag"].tolist()

    df = pd.DataFrame(index=locus_tags)
    df["Product"] = cds_meta.set_index("Locus Tag")["Product"].reindex(locus_tags).fillna("NA")
    df["Gene"] = cds_meta.set_index("Locus Tag")["Gene"].reindex(locus_tags).fillna("NA")

    # Species info
    mag_row = MAG_info[MAG_info["mag"] == MAG]
    if not mag_row.empty:
        species_name = mag_row["species_name"].values[0]
        org_id = str(mag_row["organism_id"].values[0])
        parts = species_name.split(" ")
        genus = parts[0]
        species = parts[1] if len(parts) > 1 else ""
    else:
        genus = species = species_name = org_id = "Unknown"

    genus_abbr = genus[:3].upper() if genus else "UNK"
    species_abbr = species[:2].upper() if species else "UNK"

    df["Species_Locus_Tag"] = [f"{lt}_{genus_abbr}{species_abbr}" for lt in locus_tags]
    df["Species_Name"] = species_name
    df["Organism_Identifier"] = org_id

    # GO & KEGG extraction
    go_terms = []
    kegg_terms = []

    if "DbXrefs" in cds_meta.columns:
        for x in cds_meta["DbXrefs"]:
            gos, keggs = extract_go_kegg(x)
            go_terms.append(gos)
            kegg_terms.append(keggs)
    else:
        go_terms = [[] for _ in locus_tags]
        kegg_terms = [[] for _ in locus_tags]

    df["GO_Terms"] = go_terms
    df["KEGG_Terms"] = kegg_terms

    return df



def replace_faa_headers(MAG_dir, MAG, df):
    """
    Replace FASTA headers with UniProt-style annotation.
    """
    if df is None:
        return None

    faa_in = os.path.join(MAG_dir, f"{MAG}_cds.faa")
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
                # Even though the database is not actually tr, DIA-NN can't parse it if bakta is added there
                # Tricking it into working using tr for simplicity, even though it is not strictly speaking correct. 
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



def concatenate_uniprot_fastas(base_dir, output_path):
    """
    Combines all MAG uniprot fastas into one file.
    """
    files = []

    for MAG in os.listdir(base_dir):
        f = os.path.join(base_dir, MAG, f"{MAG}_uniprot.faa")
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
    """
    Write GO and KEGG annotations to Snakemake-designated output paths.
    """
    df[["GO_Terms"]].to_csv(GO_OUT, sep="\t", header=False)
    df[["KEGG_Terms"]].to_csv(KEGG_OUT, sep="\t", header=False)

    logprint(f"Wrote GO annotations → {GO_OUT}")
    logprint(f"Wrote KEGG annotations → {KEGG_OUT}")



# -----------------------------
# MAIN
# -----------------------------
def main():
    MAG_info = pd.read_csv(MAG_METADATA_PATH, sep="\t")
    MAG_list = MAG_info["mag"].tolist()

    all_meta = []

    for MAG in MAG_list:
        MAG_dir = os.path.join(BASE_DIR, MAG)
        if not os.path.isdir(MAG_dir):
            logprint(f"[{MAG}] Skipping (missing bakta directory)")
            continue

        logprint(f"Processing {MAG}")

        preprocess_bakta_annotations(MAG_dir, MAG)
        df = format_uniprot_metadata(MAG_dir, MAG, MAG_info)
        replace_faa_headers(MAG_dir, MAG, df)

        if df is not None:
            df["MAG"] = MAG
            all_meta.append(df)

    if not all_meta:
        logprint("No MAG metadata found — nothing to export.")
        log.close()
        return

    full_df = pd.concat(all_meta)

    # Write annotation TSVs
    write_go_kegg_annotation_files(full_df)

    # Combine FASTAs
    concatenate_uniprot_fastas(BASE_DIR, FASTA_OUT)

    logprint("MAG UniProt + GO/KEGG processing complete.")
    log.close()



if __name__ == "__main__":
    main()

