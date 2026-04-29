#!/usr/bin/env python3
"""
tryptic_digest.py

Perform in-silico tryptic digestion of Prodigal-predicted proteins from all
genomes and output a peptide-to-genome mapping.

Input:
  - Prodigal .faa files (one per genome)
  - taxonomy.txt (to get genome list)

Output:
  - peptide_genome_mapping.tsv.gz with columns:
    peptide | peptide_il | genome | protein_id
"""

import gzip
import os
import re

from Bio import SeqIO
import pandas as pd

# Snakemake bindings
TAXONOMY_PATH = snakemake.input.taxonomy
PRODIGAL_DIR = snakemake.params.prodigal_dir
OUTPUT_PATH = snakemake.output.peptide_mapping
LOG_PATH = snakemake.log[0]

os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)

log = open(LOG_PATH, "w")


def logprint(msg):
    print(msg)
    print(msg, file=log, flush=True)


# Trypsin cleavage: after K or R, unless followed by P
# Returns peptides with 0 missed cleavages
TRYPSIN_PATTERN = re.compile(r"(?<=[KR])(?!P)")

MIN_PEPTIDE_LEN = 7
MAX_PEPTIDE_LEN = 30


def tryptic_digest(sequence):
    """
    In-silico tryptic digestion of a protein sequence.
    Cleaves after K/R unless followed by P. Filters by length.
    """
    # Remove stop codon marker if present
    sequence = sequence.rstrip("*")

    fragments = TRYPSIN_PATTERN.split(str(sequence))
    peptides = []
    for frag in fragments:
        if MIN_PEPTIDE_LEN <= len(frag) <= MAX_PEPTIDE_LEN:
            peptides.append(frag)
    return peptides


def main():
    # Read taxonomy to get genome list
    tax_df = pd.read_csv(TAXONOMY_PATH, sep="\t")
    genome_list = tax_df["genome"].astype(str).tolist()
    logprint(f"Found {len(genome_list)} genomes in taxonomy.txt")

    rows = []
    total_proteins = 0
    total_peptides = 0

    for genome in genome_list:
        faa_path = os.path.join(PRODIGAL_DIR, f"{genome}.faa")
        if not os.path.exists(faa_path):
            logprint(f"WARNING: No Prodigal output for genome '{genome}': {faa_path}")
            continue

        genome_proteins = 0
        genome_peptides = 0

        for record in SeqIO.parse(faa_path, "fasta"):
            genome_proteins += 1
            peptides = tryptic_digest(str(record.seq))

            for pep in peptides:
                pep_il = pep.replace("I", "L")
                rows.append((pep, pep_il, genome, record.id))
                genome_peptides += 1

        total_proteins += genome_proteins
        total_peptides += genome_peptides
        logprint(
            f"  {genome}: {genome_proteins} proteins → {genome_peptides} tryptic peptides"
        )

    logprint(f"Total: {total_proteins} proteins, {total_peptides} peptides")

    if not rows:
        logprint("WARNING: No peptides generated from any genome")
        # Write empty file with header
        with gzip.open(OUTPUT_PATH, "wt") as f:
            f.write("peptide\tpeptide_il\tgenome\tprotein_id\n")
        log.close()
        return

    df = pd.DataFrame(rows, columns=["peptide", "peptide_il", "genome", "protein_id"])
    logprint(f"Writing {len(df)} rows to {OUTPUT_PATH}")
    df.to_csv(OUTPUT_PATH, sep="\t", index=False, compression="gzip")

    logprint("Tryptic digest complete")
    log.close()


main()
