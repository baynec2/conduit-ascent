#!/usr/bin/env python3
"""
compute_peptide_lca.py

Compute the Lowest Common Ancestor (LCA) for each tryptic peptide based on the
user-provided taxonomy of the genomes that produce it. Then generate rank-specific
TSVs and FASTAs matching the peptidotyping format.

Input:
  - peptide_genome_mapping.tsv.gz (from tryptic_digest.py)
  - taxonomy.txt (user-provided lineage for each genome)

Output:
  - {rank}_lca_filtered_peptides.tsv  (×3: family, genus, species_strain)
  - {rank}_peptidotyping_db.fasta     (×3)
  - taxid_to_family_genus.tsv         (name-based lineage mapping)
"""

import csv
import os
import shlex
import subprocess

# Snakemake bindings
PEPTIDE_MAPPING_PATH = snakemake.input.peptide_mapping
TAXONOMY_PATH = snakemake.input.taxonomy

FAMILY_TSV = snakemake.output.family_tsv
FAMILY_FASTA = snakemake.output.family_fasta
GENUS_TSV = snakemake.output.genus_tsv
GENUS_FASTA = snakemake.output.genus_fasta
SPECIES_TSV = snakemake.output.species_tsv
SPECIES_FASTA = snakemake.output.species_fasta
TAXID_FAMILY_MAP = snakemake.output.taxid_family_map

LOG_PATH = snakemake.log[0]

os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
log = open(LOG_PATH, "w")


def logprint(msg):
    print(msg)
    print(msg, file=log, flush=True)


# Taxonomic ranks in order from broadest to most specific
RANKS = ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"]

# Ranks we produce databases for, and what taxonomy ranks they include
RANK_FILTERS = {
    "family": {"family"},
    "genus": {"genus"},
    "species_strain": {"species"},
}


def underscore_name(name):
    """Replace spaces with underscores and strip for use in FASTA headers."""
    if name is None:
        return "NA"
    s = str(name).strip()
    if s == "" or s.upper() == "NA":
        return "NA"
    return s.replace(" ", "_")


def compute_lca(lineages):
    """
    Compute LCA from a list of lineage dicts.
    Walk from domain → species, stop at the first rank where lineages disagree.
    Returns (lca_rank, lca_name, parent_rank, parent_name).
    """
    lca_rank = None
    lca_name = None
    parent_rank = None
    parent_name = None

    for rank in RANKS:
        values = set()
        for lin in lineages:
            val = str(lin.get(rank, "NA")).strip()
            if val.upper() == "NA" or val == "":
                val = "NA"
            values.add(val)

        # Remove NA — if all are NA, skip this rank
        values.discard("NA")

        if len(values) == 0:
            # All NA at this rank — skip, keep previous LCA
            continue
        elif len(values) == 1:
            # All agree at this rank — update LCA
            parent_rank = lca_rank
            parent_name = lca_name
            lca_rank = rank
            lca_name = values.pop()
        else:
            # Disagreement — stop here
            break

    return lca_rank, lca_name, parent_rank, parent_name


def emit_peptide(pep_il, rep_peptide, genomes, genome_lineage, tsv_files,
                 fasta_files, counts, pid):
    """Compute the LCA for one peptide's genome set and write its rank outputs.

    Returns True if the peptide consumed the peptide id `pid` (i.e. it had a
    valid LCA), False if it was skipped for lack of lineage. Whether a row is
    actually written depends on the LCA rank (only family/genus/species ranks
    are useful for detection); higher-level LCAs still consume an id, matching
    the original numbering.
    """
    lineages = [genome_lineage[g] for g in genomes if g in genome_lineage]
    if not lineages:
        return False

    lca_rank, lca_name, parent_rank, parent_name = compute_lca(lineages)
    if lca_rank is None:
        # No agreement at any rank (all NA) — skip
        return False

    lca_name_u = underscore_name(lca_name)
    parent_name_u = underscore_name(parent_name if parent_name else "NA")

    # Map the LCA rank to our database rank categories
    if lca_rank == "family":
        db_rank = "family"
    elif lca_rank == "genus":
        db_rank = "genus"
    elif lca_rank in ("species", "strain"):
        db_rank = "species_strain"
    else:
        # Higher-level LCA (phylum, class, order, etc.) — id consumed, no row
        return True

    fasta_header = (
        f"gpep|{pid}|{lca_name_u} {lca_rank}_{lca_name_u} "
        f"OS={lca_name} OX={lca_name_u} RK={lca_rank} PT={parent_name_u}"
    )

    # TSV row (10 columns matching peptidotyping format); lca/lca_il and fa/fa_il
    # are the same name-based identifier.
    tsv_files[db_rank].write(
        f"{pid}\t{rep_peptide}\t{lca_name_u}\t{lca_name_u}\t"
        f"{lca_name_u}\t{lca_name_u}\t{lca_name}\t{lca_rank}\t"
        f"{parent_name_u}\t{fasta_header}\n"
    )
    fasta_files[db_rank].write(f">{fasta_header}\n{rep_peptide}\n")
    counts[db_rank] += 1
    return True


def main():
    logprint("Loading taxonomy")
    # Stream taxonomy.txt (one row per genome) into a genome → lineage dict.
    # This is small (bounded by genome count); no need for pandas.
    genome_lineage = {}
    with open(TAXONOMY_PATH, newline="") as tf:
        reader = csv.DictReader(tf, delimiter="\t")
        for row in reader:
            genome = str(row["genome"])
            lineage = {}
            for rank in RANKS:
                val = row.get(rank)
                lineage[rank] = str(val).strip() if val is not None else "NA"
            genome_lineage[genome] = lineage
    logprint(f"  {len(genome_lineage)} genomes in taxonomy")

    # ── Build taxid_to_family_genus.tsv ──────────────────────────────────────
    # For every unique taxon that appears in a lineage, record its family and
    # genus. Used by infer_family_presence.R to map lca → family. Depends only
    # on the taxonomy, not the peptide mapping.
    logprint("Building taxid_to_family_genus.tsv")
    all_taxa = {}  # name_u → (rank, family_u, genus_u)
    for genome, lineage in genome_lineage.items():
        family_u = underscore_name(lineage.get("family", "NA"))
        genus_u = underscore_name(lineage.get("genus", "NA"))
        for rank in RANKS:
            name_u = underscore_name(lineage.get(rank, "NA"))
            if name_u == "NA":
                continue
            if name_u not in all_taxa:
                all_taxa[name_u] = (rank, family_u, genus_u)

    with open(TAXID_FAMILY_MAP, "w") as f:
        for name_u, (rank, family_u, genus_u) in sorted(all_taxa.items()):
            f.write(f"{name_u}\t{rank}\t{family_u}\t{genus_u}\n")
    logprint(f"  {len(all_taxa)} entries in taxid_to_family_genus.tsv")

    # ── Open rank-specific outputs ───────────────────────────────────────────
    TSV_HEADER = (
        "id\tsequence\tlca\tlca_il\tfa\tfa_il\tname\trank\tparent_id\tfasta_header\n"
    )
    output_map = {
        "family": (FAMILY_TSV, FAMILY_FASTA),
        "genus": (GENUS_TSV, GENUS_FASTA),
        "species_strain": (SPECIES_TSV, SPECIES_FASTA),
    }
    counts = {"family": 0, "genus": 0, "species_strain": 0}
    tsv_files = {}
    fasta_files = {}
    for db_rank, (tsv_path, fasta_path) in output_map.items():
        os.makedirs(os.path.dirname(tsv_path), exist_ok=True)
        tsv_files[db_rank] = open(tsv_path, "w")
        tsv_files[db_rank].write(TSV_HEADER)
        fasta_files[db_rank] = open(fasta_path, "w")

    # ── Stream peptides grouped by peptide_il ────────────────────────────────
    # The mapping is far too large to hold in memory (a full catalog is
    # 10s-of-millions of rows), so we group by peptide via an external disk-based
    # `sort` on the peptide_il column (field 2) and walk the sorted stream one
    # peptide at a time. Memory stays bounded by the genome set of a single
    # peptide. `sort` uses $TMPDIR (Snakemake points it at cluster scratch).
    logprint("Sorting peptide mapping by peptide and computing LCA")
    tmpdir = os.environ.get("TMPDIR", os.path.dirname(FAMILY_TSV))
    sort_pipeline = (
        "set -o pipefail; "
        f"gzip -dc {shlex.quote(PEPTIDE_MAPPING_PATH)} | tail -n +2 | "
        f"LC_ALL=C sort -t '\t' -k2,2 -T {shlex.quote(tmpdir)}"
    )
    proc = subprocess.Popen(
        sort_pipeline, shell=True, executable="/bin/bash",
        stdout=subprocess.PIPE, text=True,
    )

    peptide_id = 0
    skipped_no_lineage = 0
    n_peptides = 0
    cur_pep_il = None
    rep_peptide = None
    genomes = set()

    for line in proc.stdout:
        line = line.rstrip("\n")
        if not line:
            continue
        # peptide, peptide_il, genome, protein_id
        peptide, pep_il, genome = line.split("\t")[:3]
        if pep_il != cur_pep_il:
            if cur_pep_il is not None:
                n_peptides += 1
                pid = peptide_id + 1
                if emit_peptide(cur_pep_il, rep_peptide, genomes, genome_lineage,
                                tsv_files, fasta_files, counts, pid):
                    peptide_id = pid
                else:
                    skipped_no_lineage += 1
            cur_pep_il = pep_il
            rep_peptide = peptide  # first (sorted) original sequence for this peptide
            genomes = set()
        genomes.add(genome)

    # Flush the final peptide group
    if cur_pep_il is not None:
        n_peptides += 1
        pid = peptide_id + 1
        if emit_peptide(cur_pep_il, rep_peptide, genomes, genome_lineage,
                        tsv_files, fasta_files, counts, pid):
            peptide_id = pid
        else:
            skipped_no_lineage += 1

    proc.stdout.close()
    if proc.wait() != 0:
        for fh in list(tsv_files.values()) + list(fasta_files.values()):
            fh.close()
        raise RuntimeError(f"sort pipeline failed with exit code {proc.returncode}")

    logprint(
        f"  {n_peptides} unique peptides; {peptide_id} with LCA assigned "
        f"({skipped_no_lineage} skipped for no lineage)"
    )

    for db_rank in output_map:
        tsv_files[db_rank].close()
        fasta_files[db_rank].close()

    logprint("Rank database generation complete:")
    for db_rank, count in counts.items():
        logprint(f"  {db_rank}: {count} peptides")

    logprint("compute_peptide_lca.py finished")
    log.close()


main()
