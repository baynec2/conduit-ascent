#!/usr/bin/env python3
"""
build_marker_gene_fasta.py

Parse hmmscan --tblout output files for each genome and extract the
corresponding marker gene (ribP/elonF) protein sequences from FragGeneScan
FAA files. Write all hits to a combined FASTA with headers:
  >{genome_id}|{protein_id}

Premature stop codons (*) are stripped from sequences.
"""

import os
import sys
from Bio import SeqIO

# Snakemake bindings
HMMER_FILES = snakemake.input.hmmer_files   # list of .txt files, one per genome
FAA_FILES   = snakemake.input.faa_files     # list of .faa files, one per genome
OUT_FASTA   = snakemake.output[0]
LOG         = snakemake.log[0]

os.makedirs(os.path.dirname(OUT_FASTA), exist_ok=True)

log = open(LOG, "w")
def logprint(msg):
    print(msg)
    print(msg, file=log)


def parse_hmmscan_tblout(hmmscan_file, evalue_threshold=1e-10):
    """
    Parse hmmscan --tblout output. Return set of query protein IDs with
    at least one hit below evalue_threshold.
    Column indices (0-based): 0=target, 2=query, 4=full_evalue
    """
    hits = set()
    with open(hmmscan_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) < 5:
                continue
            try:
                evalue = float(cols[4])
            except ValueError:
                continue
            if evalue <= evalue_threshold:
                hits.add(cols[2])  # query = protein ID
    return hits


def genome_id_from_hmmer_path(hmmer_path):
    """Extract genome ID from path like hapid/hmmer/{genome}_hmmer.txt"""
    basename = os.path.basename(hmmer_path)
    return basename.replace("_hmmer.txt", "")


# Build a map from genome_id → faa_path
faa_map = {}
for faa_path in FAA_FILES:
    gid = os.path.splitext(os.path.basename(faa_path))[0]
    faa_map[gid] = faa_path

total_hits = 0
with open(OUT_FASTA, "w") as out_fh:
    for hmmer_file in HMMER_FILES:
        genome_id = genome_id_from_hmmer_path(hmmer_file)
        hit_ids = parse_hmmscan_tblout(hmmer_file)
        logprint(f"[{genome_id}] {len(hit_ids)} marker gene hits")

        faa_path = faa_map.get(genome_id)
        if faa_path is None:
            logprint(f"[{genome_id}] WARNING: no FAA file found — skipping")
            continue

        for record in SeqIO.parse(faa_path, "fasta"):
            if record.id in hit_ids:
                seq_str = str(record.seq).rstrip("*")  # strip premature stop
                out_fh.write(f">{genome_id}|{record.id}\n{seq_str}\n")
                total_hits += 1

logprint(f"Total marker gene sequences written: {total_hits}")
log.close()
