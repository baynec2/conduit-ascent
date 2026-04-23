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

import os

import pandas as pd

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
    if pd.isna(name) or str(name).strip().upper() == "NA":
        return "NA"
    return str(name).strip().replace(" ", "_")


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


def main():
    logprint("Loading taxonomy")
    tax_df = pd.read_csv(TAXONOMY_PATH, sep="\t")
    logprint(f"  {len(tax_df)} genomes in taxonomy")

    # Build genome → lineage dict
    genome_lineage = {}
    for _, row in tax_df.iterrows():
        genome = str(row["genome"])
        lineage = {}
        for rank in RANKS:
            if rank in row:
                lineage[rank] = str(row[rank]).strip()
            else:
                lineage[rank] = "NA"
        genome_lineage[genome] = lineage

    logprint("Loading peptide-genome mapping")
    pep_df = pd.read_csv(PEPTIDE_MAPPING_PATH, sep="\t", compression="gzip")
    logprint(f"  {len(pep_df)} rows in peptide mapping")

    # Group by peptide_il to find all genomes producing each peptide
    logprint("Grouping peptides and computing LCA")
    grouped = pep_df.groupby("peptide_il")

    peptide_records = []
    peptide_id = 0
    skipped_no_lineage = 0

    for pep_il, group in grouped:
        genomes = group["genome"].unique().tolist()
        # Pick one representative original peptide sequence
        rep_peptide = group["peptide"].iloc[0]

        # Get lineages for all genomes
        lineages = []
        for g in genomes:
            if g in genome_lineage:
                lineages.append(genome_lineage[g])
        if not lineages:
            skipped_no_lineage += 1
            continue

        lca_rank, lca_name, parent_rank, parent_name = compute_lca(lineages)
        if lca_rank is None:
            # No agreement at any rank (all NA) — skip
            skipped_no_lineage += 1
            continue

        peptide_id += 1
        peptide_records.append(
            {
                "id": peptide_id,
                "peptide": rep_peptide,
                "peptide_il": pep_il,
                "lca_rank": lca_rank,
                "lca_name": lca_name,
                "parent_rank": parent_rank,
                "parent_name": parent_name,
            }
        )

    logprint(
        f"  {len(peptide_records)} unique peptides with LCA assigned "
        f"({skipped_no_lineage} skipped)"
    )

    # ── Build taxid_to_family_genus.tsv ──────────────────────────────────────
    # For every unique taxon that appears as an LCA, record its family and genus.
    # This is used by infer_family_presence.R to map lca → family.
    logprint("Building taxid_to_family_genus.tsv")

    # Collect all unique taxa from lineages
    all_taxa = {}  # name_u → (rank, family_u, genus_u)
    for genome, lineage in genome_lineage.items():
        family_u = underscore_name(lineage.get("family", "NA"))
        genus_u = underscore_name(lineage.get("genus", "NA"))
        for rank in RANKS:
            name = lineage.get(rank, "NA")
            name_u = underscore_name(name)
            if name_u == "NA":
                continue
            if name_u not in all_taxa:
                all_taxa[name_u] = (rank, family_u, genus_u)

    with open(TAXID_FAMILY_MAP, "w") as f:
        for name_u, (rank, family_u, genus_u) in sorted(all_taxa.items()):
            f.write(f"{name_u}\t{rank}\t{family_u}\t{genus_u}\n")

    logprint(f"  {len(all_taxa)} entries in taxid_to_family_genus.tsv")

    # ── Generate rank-specific TSVs and FASTAs ───────────────────────────────
    TSV_HEADER = (
        "id\tsequence\tlca\tlca_il\tfa\tfa_il\tname\trank\tparent_id\tfasta_header\n"
    )

    output_map = {
        "family": (FAMILY_TSV, FAMILY_FASTA),
        "genus": (GENUS_TSV, GENUS_FASTA),
        "species_strain": (SPECIES_TSV, SPECIES_FASTA),
    }

    counts = {"family": 0, "genus": 0, "species_strain": 0}

    # Open all output files
    tsv_files = {}
    fasta_files = {}
    for db_rank, (tsv_path, fasta_path) in output_map.items():
        os.makedirs(os.path.dirname(tsv_path), exist_ok=True)
        tsv_files[db_rank] = open(tsv_path, "w")
        tsv_files[db_rank].write(TSV_HEADER)
        fasta_files[db_rank] = open(fasta_path, "w")

    for rec in peptide_records:
        lca_rank = rec["lca_rank"]
        lca_name = rec["lca_name"]
        lca_name_u = underscore_name(lca_name)
        parent_name_u = underscore_name(
            rec["parent_name"] if rec["parent_name"] else "NA"
        )

        # Map the LCA rank to our database rank categories
        if lca_rank == "family":
            db_rank = "family"
        elif lca_rank == "genus":
            db_rank = "genus"
        elif lca_rank in ("species", "strain"):
            db_rank = "species_strain"
        else:
            # Higher-level LCA (phylum, class, order, etc.) — not useful for detection
            continue

        pid = rec["id"]
        sequence = rec["peptide"]

        # Build FASTA header matching peptidotyping format
        fasta_header = (
            f"gp|{pid}|{lca_name_u} {lca_rank}_{lca_name_u} "
            f"OS={lca_name} OX={lca_name_u} RK={lca_rank} PT={parent_name_u}"
        )

        # Write TSV row (10 columns matching peptidotyping format)
        # lca and lca_il are the same (name-based identifier)
        # fa and fa_il are set to the same value
        tsv_files[db_rank].write(
            f"{pid}\t{sequence}\t{lca_name_u}\t{lca_name_u}\t"
            f"{lca_name_u}\t{lca_name_u}\t{lca_name}\t{lca_rank}\t"
            f"{parent_name_u}\t{fasta_header}\n"
        )

        # Write FASTA
        fasta_files[db_rank].write(f">{fasta_header}\n{sequence}\n")

        counts[db_rank] += 1

    for db_rank in output_map:
        tsv_files[db_rank].close()
        fasta_files[db_rank].close()

    logprint("Rank database generation complete:")
    for db_rank, count in counts.items():
        logprint(f"  {db_rank}: {count} peptides")

    logprint("compute_peptide_lca.py finished")
    log.close()


main()
