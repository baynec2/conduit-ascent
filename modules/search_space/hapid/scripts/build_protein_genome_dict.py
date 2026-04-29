#!/usr/bin/env python3
"""
build_protein_genome_dict.py

Parse the deduplicated marker gene FASTA (after CD-HIT) and build a
protein→genome JSON dictionary from the {genome_id}|{protein_id} headers.

Output: protein2genome_dic.json  { protein_id: genome_id }
"""

import json
import os
from Bio import SeqIO

# Snakemake bindings
FASTA_IN = snakemake.input.fasta
JSON_OUT  = snakemake.output[0]
LOG       = snakemake.log[0]

os.makedirs(os.path.dirname(JSON_OUT), exist_ok=True)

log = open(LOG, "w")
def logprint(msg):
    print(msg)
    print(msg, file=log)

protein2genome = {}
for record in SeqIO.parse(FASTA_IN, "fasta"):
    header = record.id  # genome_id|protein_id
    if "|" not in header:
        logprint(f"WARNING: unexpected header format (no pipe): {header}")
        continue
    genome_id, protein_id = header.split("|", 1)
    protein2genome[protein_id] = genome_id

logprint(f"Built protein→genome dict with {len(protein2genome)} entries")

with open(JSON_OUT, "w") as f:
    json.dump(protein2genome, f)

logprint(f"Written to {JSON_OUT}")
log.close()
