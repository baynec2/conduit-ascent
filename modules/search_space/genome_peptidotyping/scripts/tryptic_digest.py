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
import csv

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
    # Read taxonomy to get genome list (first column is "genome"). Stream the
    # file instead of loading it into pandas — we only need the genome column.
    genome_list = []
    with open(TAXONOMY_PATH, newline="") as tax:
        reader = csv.DictReader(tax, delimiter="\t")
        for rec in reader:
            genome_list.append(str(rec["genome"]))
    logprint(f"Found {len(genome_list)} genomes in taxonomy.txt")

    total_proteins = 0
    total_peptides = 0

    # Stream peptides straight to the gzip output as they are produced. Memory
    # stays flat (one genome at a time) regardless of catalog size — the old
    # approach accumulated every peptide of every genome into a list + pandas
    # DataFrame, which OOM'd on full catalogs (thousands of genomes → tens of
    # millions of peptides).
    with gzip.open(OUTPUT_PATH, "wt", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["peptide", "peptide_il", "genome", "protein_id"])

        for genome in genome_list:
            faa_path = os.path.join(PRODIGAL_DIR, f"{genome}.faa")
            if not os.path.exists(faa_path):
                logprint(f"WARNING: No Prodigal output for genome '{genome}': {faa_path}")
                continue

            genome_proteins = 0
            genome_peptides = 0

            for record in SeqIO.parse(faa_path, "fasta"):
                genome_proteins += 1
                for pep in tryptic_digest(str(record.seq)):
                    writer.writerow((pep, pep.replace("I", "L"), genome, record.id))
                    genome_peptides += 1

            total_proteins += genome_proteins
            total_peptides += genome_peptides
            logprint(
                f"  {genome}: {genome_proteins} proteins → {genome_peptides} tryptic peptides"
            )

    logprint(f"Total: {total_proteins} proteins, {total_peptides} peptides")
    if total_peptides == 0:
        logprint("WARNING: No peptides generated from any genome")
    logprint("Tryptic digest complete")
    log.close()


main()
