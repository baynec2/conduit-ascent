#!/usr/bin/env bash
# Tests the pure-awk transforms inside modules/search_space/unipept_peptidotyping/
# unipept_peptidotyping.smk. The rule-internal awk scripts are reproduced here
# verbatim except for snakemake's double-brace escaping; the .smk uses {{...}}
# inside shell directives where this file uses {...}. If you change the awk in
# the .smk, mirror the change here (and vice versa) — there is no shared source.
#
# Covered:
#   1. effective-rank-priority assignment (family > genus > species_strain,
#      finest rank with aggregate peptide count >= min_peptides)
#   2. FAM= tag injection on family-rank FASTA entries
#   3. FAM= tag injection on genus-rank FASTA entries via the lineage join
#
# Usage: bash tests/scripts/test_awk_pipelines.sh
# Exits 0 on success, 1 on failure.

set -euo pipefail

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

PASS=0
FAIL=0

# Pretty-print test outcome and bump counters.
report() {
    local name="$1" status="$2"
    if [ "$status" = "pass" ]; then
        printf "  PASS  %s\n" "$name"
        PASS=$((PASS + 1))
    else
        printf "  FAIL  %s\n" "$name"
        FAIL=$((FAIL + 1))
    fi
}

# Diff `actual` against `expected`; report by name.
assert_eq_file() {
    local name="$1" expected="$2" actual="$3"
    if diff -u "$expected" "$actual" > "$TMP/$name.diff"; then
        report "$name" pass
    else
        report "$name" fail
        echo "    --- diff (expected vs actual) ---"
        sed 's/^/    /' "$TMP/$name.diff"
    fi
}

# ─────────────────────────────────────────────────────────────────────────────
# Test 1: rank-priority assignment
#
# Input COUNTS_FILE rows (family_taxid, rank, n_peptides). With min=10:
#   F1 has family:15, genus:50, species_strain:200 → finest passing = family
#   F2 has family:5,  genus:50, species_strain:200 → finest passing = genus
#   F3 has family:5,  genus:5,  species_strain:200 → finest passing = species_strain
#   F4 has family:5,  genus:5,  species_strain:5   → excluded (none pass min)
# ─────────────────────────────────────────────────────────────────────────────

cat > "$TMP/counts.tsv" <<'EOF'
family_taxid	rank	n_peptides
F1	family	15
F1	genus	50
F1	species_strain	200
F2	family	5
F2	genus	50
F2	species_strain	200
F3	family	5
F3	genus	5
F3	species_strain	200
F4	family	5
F4	genus	5
F4	species_strain	5
EOF

# Mirrors unipept_peptidotyping.smk lines 327-349 (single-brace form).
awk -F'\t' -v min="10" '
    NR==1 { next }
    {
        fam=$1; rank=$2; n=$3
        count[fam"\t"rank] = n
        families[fam] = 1
    }
    END {
        for (fam in families) {
            fam_key = fam"\tfamily"
            gen_key = fam"\tgenus"
            sps_key = fam"\tspecies_strain"
            if (fam_key in count && count[fam_key]+0 >= min+0) {
                print fam"\tfamily\t"count[fam_key]
            } else if (gen_key in count && count[gen_key]+0 >= min+0) {
                print fam"\tgenus\t"count[gen_key]
            } else if (sps_key in count && count[sps_key]+0 >= min+0) {
                print fam"\tspecies_strain\t"count[sps_key]
            }
        }
    }
' "$TMP/counts.tsv" | sort > "$TMP/rank_mapping.actual"

# Sort the expected output since awk iterates families in arbitrary order.
cat > "$TMP/rank_mapping.expected" <<'EOF'
F1	family	15
F2	genus	50
F3	species_strain	200
EOF
sort -o "$TMP/rank_mapping.expected" "$TMP/rank_mapping.expected"

assert_eq_file "rank-priority assignment" "$TMP/rank_mapping.expected" "$TMP/rank_mapping.actual"

# ─────────────────────────────────────────────────────────────────────────────
# Test 2: family-rank FASTA emission with FAM= tag
#
# family_tsv rows have columns:
#   id sequence lca lca_il fa fa_il name rank parent_id fasta_header
# When the family (lca_il) appears in rank_lookup with effective_rank==family,
# emit ">$fasta_header FAM=$lca_il" followed by the sequence.
# ─────────────────────────────────────────────────────────────────────────────

cat > "$TMP/family_tsv" <<'EOF'
id	sequence	lca	lca_il	fa	fa_il	name	rank	parent_id	fasta_header
1	AAAAK	Enterobacteriaceae	543	-	-	Enterobacteriaceae	family	91347	umgap|1|543 family_Enterobacteriaceae OS=Enterobacteriaceae OX=543 RK=family PT=91347
2	BBBBK	Bacillaceae	186817	-	-	Bacillaceae	family	1385	umgap|2|186817 family_Bacillaceae OS=Bacillaceae OX=186817 RK=family PT=1385
EOF

cat > "$TMP/rank_lookup" <<'EOF'
543	family
186817	genus
EOF

# Mirrors unipept_peptidotyping.smk lines 368-381 (single-brace form).
awk -F'\t' '
    NR==FNR {
        if ($2=="family") fam_families[$1]=1
        next
    }
    FNR==1 { next }
    $4 in fam_families {
        fam_taxid = $4
        header = $10 " FAM=" fam_taxid
        print ">" header
        print $2
    }
' "$TMP/rank_lookup" "$TMP/family_tsv" > "$TMP/family_fasta.actual"

cat > "$TMP/family_fasta.expected" <<'EOF'
>umgap|1|543 family_Enterobacteriaceae OS=Enterobacteriaceae OX=543 RK=family PT=91347 FAM=543
AAAAK
EOF

assert_eq_file "family-rank FASTA + FAM=" "$TMP/family_fasta.expected" "$TMP/family_fasta.actual"

# ─────────────────────────────────────────────────────────────────────────────
# Test 3: genus-rank FASTA emission via lineage join
#
# genus_tsv lca_il is a genus taxid. We join against the lineage file to
# recover the parent family_taxid, then require that family to be assigned
# the "genus" effective rank. Output gets the parent family as FAM=.
#
# Three-file awk: lineage → rank_lookup → genus_tsv.
# ─────────────────────────────────────────────────────────────────────────────

cat > "$TMP/lineage" <<'EOF'
561	genus	543	561
1386	genus	186817	1386
562	species	543	561
EOF

cat > "$TMP/rank_lookup2" <<'EOF'
543	family
186817	genus
EOF

cat > "$TMP/genus_tsv" <<'EOF'
id	sequence	lca	lca_il	fa	fa_il	name	rank	parent_id	fasta_header
1	GGGGK	Escherichia	561	-	-	Escherichia	genus	543	umgap|1|561 genus_Escherichia OS=Escherichia OX=561 RK=genus PT=543
2	HHHHK	Bacillus	1386	-	-	Bacillus	genus	186817	umgap|2|1386 genus_Bacillus OS=Bacillus OX=1386 RK=genus PT=186817
EOF

# Mirrors unipept_peptidotyping.smk lines 401-417 (single-brace form).
awk -F'\t' '
    ARGIND==1 { fam[$1]=$3; next }
    ARGIND==2 { if($2=="genus") genus_fams[$1]=1; next }
    FNR==1 { next }
    {
        lca_il=$4
        if (lca_il in fam && fam[lca_il] in genus_fams) {
            fam_taxid = fam[lca_il]
            header = $10 " FAM=" fam_taxid
            print ">" header
            print $2
        }
    }
' "$TMP/lineage" "$TMP/rank_lookup2" "$TMP/genus_tsv" > "$TMP/genus_fasta.actual"

# Only the Bacillus (1386) entry should be emitted:
#   - 1386's parent family is 186817, which is assigned "genus" rank → kept
#   - 561's parent family is 543, which is assigned "family" rank → excluded
cat > "$TMP/genus_fasta.expected" <<'EOF'
>umgap|2|1386 genus_Bacillus OS=Bacillus OX=1386 RK=genus PT=186817 FAM=186817
HHHHK
EOF

assert_eq_file "genus-rank FASTA + FAM= via lineage join" "$TMP/genus_fasta.expected" "$TMP/genus_fasta.actual"

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────

echo ""
echo "Results: $PASS passed, $FAIL failed"
[ "$FAIL" -eq 0 ]
