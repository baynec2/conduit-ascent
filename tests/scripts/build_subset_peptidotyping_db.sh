#!/usr/bin/env bash
# Build a tiny subset of the peptidotyping per-rank peptide TSVs for smoke
# tests. The subset is filtered to a curated allowlist of taxa, capped at
# MAX_PER_TAXON peptides per taxon at each rank, and written to
# tests/fixtures/peptidotyping_subset/ so a smoke run can point
# config["peptidotyping_resource_dir"] there.
#
# Inputs (from the full resources, default location is the host profile's
# peptidotyping_resource_dir):
#   family_lca_filtered_peptides.tsv
#   genus_lca_filtered_peptides.tsv
#   species_strain_lca_filtered_peptides.tsv
#   taxons.tsv.lz4    (passed through, no subsetting needed — already small)
#
# Output:
#   tests/fixtures/peptidotyping_subset/{family,genus,species_strain}_lca_filtered_peptides.tsv
#   tests/fixtures/peptidotyping_subset/allowlist.txt   (the curated taxid list, checked in)
#
# Usage:
#   bash tests/scripts/build_subset_peptidotyping_db.sh \
#        /home/.../HDD/peptidotyping_resources/   # SOURCE_DIR (defaults to host)
#
# The species_strain TSV is large (~392 GB on disk). To keep the build
# tractable, we early-exit once every allowlisted taxon has hit its cap; this
# typically means we never need to scan past the first few GB.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SOURCE_DIR="${1:-/home/nanopore-catalyst/HDD/peptidotyping_resources}"
FIXTURE_DIR="$REPO_ROOT/tests/fixtures/peptidotyping_subset"

# How many peptides to keep per taxon at each rank. Smoke peptidotyping needs
# enough to clear FDR for the target organism; ~1000-2000 is plenty.
MAX_PER_TAXON="${MAX_PER_TAXON:-2000}"

mkdir -p "$FIXTURE_DIR"

# ── Curated allowlist ─────────────────────────────────────────────────────────
# The 10 species in the ATCC-defined-community pool used for benchmarking
# (see experiments/pool_pilot/input/ncbi_taxa_ids.txt) plus the genus and
# family taxid for each — so the rank-priority logic in
# build_effective_detection_rank_db has data at every level.
#
# Format: one taxid per line. Comments allowed.

cat > "$FIXTURE_DIR/allowlist.txt" <<'EOF'
# ── Pool species (10) ──────────────────────────────────────────────────
820       # Bacteroides uniformis           (genus 816, family 815)
562       # Escherichia coli                (genus 561, family 543)
1579      # Lactobacillus acidophilus       (genus 1578, family 33958)
1624      # Ligilactobacillus salivarius    (genus 2767887, family 33958)
33035     # Blautia producta                (genus 572511, family 186803)
823       # Parabacteroides distasonis      (genus 375288, family 2005525)
1679      # Bifidobacterium longum          (genus 1678, family 31953)
239935    # Akkermansia muciniphila         (genus 239934, family 1647988)
821       # Phocaeicola vulgatus            (genus 909656, family 815)
817       # Bacteroides fragilis            (genus 816, family 815)

# ── Pool genera (9 unique) ─────────────────────────────────────────────
816       # Bacteroides
561       # Escherichia
1578      # Lactobacillus
2767887   # Ligilactobacillus
572511    # Blautia
375288    # Parabacteroides
1678      # Bifidobacterium
239934    # Akkermansia
909656    # Phocaeicola

# ── Pool families (7 unique) ───────────────────────────────────────────
815       # Bacteroidaceae        (covers B. uniformis, P. vulgatus, B. fragilis)
543       # Enterobacteriaceae    (covers E. coli)
33958     # Lactobacillaceae      (covers L. acidophilus, L. salivarius)
186803    # Lachnospiraceae       (covers Blautia producta)
2005525   # Tannerellaceae        (covers Parabacteroides distasonis)
31953     # Bifidobacteriaceae    (covers Bifidobacterium longum)
1647988   # Akkermansiaceae       (covers Akkermansia muciniphila)
EOF

# Strip comments and take only the first whitespace-separated token (the
# taxid). Bare `sed 's/#.*//'` would leave trailing whitespace on every line,
# which would then mis-key the awk filter (allowed["820   "] != allowed["820"]).
ALLOW_TXT="$FIXTURE_DIR/.allowlist.plain"
sed 's/#.*//' "$FIXTURE_DIR/allowlist.txt" | awk 'NF { print $1 }' > "$ALLOW_TXT"

n_allow=$(wc -l < "$ALLOW_TXT")
echo "Allowlist contains $n_allow taxa"

# ── Subset filter (shared awk logic) ──────────────────────────────────────────
# Stream a TSV; keep rows whose col 4 (lca_il) is in the allowlist, up to
# MAX_PER_TAXON per taxon. The optional expected_taxa argument (count of
# allowlisted taxa actually expected in this rank's file) drives early-exit:
# pass 0 to scan the whole file. The species_strain TSV is ~392 GB so passing
# the expected count there saves an order of magnitude of disk I/O when the
# matching taxids cluster early.

subset_tsv() {
    local in="$1" out="$2" rank_label="$3" expected_taxa="${4:-0}"
    echo "→ filtering $rank_label rank from $(basename "$in")  (expected_taxa=$expected_taxa)"

    local t0=$(date +%s)
    awk -F'\t' -v allow="$ALLOW_TXT" -v cap="$MAX_PER_TAXON" -v expected="$expected_taxa" '
        BEGIN {
            while ((getline line < allow) > 0) {
                allowed[line] = 1
            }
        }
        NR == 1 { print; next }
        ($4 in allowed) && (kept[$4] < cap+0) {
            print
            kept[$4]++
            if (kept[$4] == cap+0) full_taxa++
            # Only early-exit when caller has told us how many taxa to expect
            # (otherwise we have no way to know we have everything).
            if (expected+0 > 0 && full_taxa == expected+0) exit
        }
    ' "$in" > "$out"

    local t1=$(date +%s)
    local nrows=$(wc -l < "$out")
    echo "  kept $((nrows - 1)) peptides  (header + data; $((t1 - t0))s)"
}

# ── Run subsetting for each rank ──────────────────────────────────────────────
# Pass the count of allowlisted taxa at each rank as the early-exit trigger.
# Allowlist has: 10 species, 9 genera (the "Pool genera" block), 7 families
# (the "Pool families" block).
subset_tsv \
    "$SOURCE_DIR/family_lca_filtered_peptides.tsv" \
    "$FIXTURE_DIR/family_lca_filtered_peptides.tsv" \
    "family" 7

subset_tsv \
    "$SOURCE_DIR/genus_lca_filtered_peptides.tsv" \
    "$FIXTURE_DIR/genus_lca_filtered_peptides.tsv" \
    "genus" 9

subset_tsv \
    "$SOURCE_DIR/species_strain_lca_filtered_peptides.tsv" \
    "$FIXTURE_DIR/species_strain_lca_filtered_peptides.tsv" \
    "species_strain" 10

# Clean up the working file; allowlist.txt remains for traceability.
rm -f "$ALLOW_TXT"

# ── Other small files the workflow consumes ──────────────────────────────────
# relnotes.txt (~1 KB): the `check_sequence_index_version` rule reads its top
# line to detect upstream UniProt version drift. Include verbatim.
if [ -f "$SOURCE_DIR/relnotes.txt" ]; then
    cp "$SOURCE_DIR/relnotes.txt" "$FIXTURE_DIR/relnotes.txt"
    echo "→ copied relnotes.txt"
fi

# taxons.tsv.lz4 (~46 MB): consumed by extract_peptidotyping_resource_metrics
# and generate_peptidotyping_db. Both rules' outputs are present in this
# fixture, so snakemake skips them — but the input declarations still resolve
# against the resource dir. We avoid copying it to keep fixture small; if a
# future rule depends on its content (not just existence), revisit.

echo
echo "Fixture written to $FIXTURE_DIR"
ls -lh "$FIXTURE_DIR"
echo
echo "Total fixture size: $(du -sh "$FIXTURE_DIR" | awk '{print $1}')"
