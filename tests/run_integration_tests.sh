#!/usr/bin/env bash
# Run end-to-end integration tests for one or all search_space_method configurations.
#
# Usage:
#   bash tests/run_integration_tests.sh [method]
#
#   method: ncbi_taxonomy_id | uniprot_proteome_id | MAGs | metaphlan | peptidotyping | all
#   Defaults to "all" if not specified.
#
# Prerequisites:
#   - snakemake >= 8
#   - apptainer/singularity
#   - tests/data/sample1.raw (tracked via Git LFS)
#
# The metaphlan method mocks the MetaPhlAn database + profiling steps via --omit-from.
# A pre-committed mock merged_profiles.txt is seeded into the run directory so that
# downstream rules (call_ncbi_taxa_ids and beyond) are not also pruned from the DAG.
#
# The MAGs method requires a FASTA file at:
#   experiments/integration_test/input/MAG_files/ecoli_mag.fa
#
# The peptidotyping method requires a pre-built sequence index at:
#   resources/peptidotyping/

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW_FILE="$REPO_ROOT/tests/data/sample1.raw"
METHOD="${1:-all}"

PASSED=()
FAILED=()

# ── Preflight checks ──────────────────────────────────────────────────────────

if [ ! -f "$RAW_FILE" ]; then
    echo "ERROR: $RAW_FILE not found."
    echo "Ensure Git LFS is installed and run: git lfs pull"
    exit 1
fi

# ── Runner ────────────────────────────────────────────────────────────────────

run_integration_test() {
    local method="$1"
    local config="$REPO_ROOT/tests/configs/integration_test_${method}.yaml"
    local extra_flags="${2:-}"

    echo "============================================================"
    echo "Integration test: $method"
    echo "============================================================"

    if snakemake \
        --snakefile "$REPO_ROOT/Snakefile" \
        --configfile "$config" \
        --use-singularity \
        --cores 4 \
        $extra_flags \
        2>&1; then
        echo "PASS: $method"
        PASSED+=("$method")
    else
        echo "FAIL: $method"
        FAILED+=("$method")
    fi
    echo
}

# ── Method dispatch ───────────────────────────────────────────────────────────

run_ncbi_taxonomy_id()    { run_integration_test ncbi_taxonomy_id; }
run_uniprot_proteome_id() {
    run_integration_test uniprot_proteome_id
    # Validate that eggnogmapper-derived annotation types are present in the
    # final conduit_annotations.txt (catches protein-ID format mismatch bugs).
    local annotations
    annotations=$(find "$REPO_ROOT/experiments/integration_test/runs/uniprot_proteome_id" \
        -name "conduit_annotations.txt" | head -1)
    if [ -n "$annotations" ]; then
        local emapper_types
        emapper_types=$(awk '{print $2}' "$annotations" | grep -v "^uniprot_" | grep -v "^annotation_type$" | sort -u | wc -l)
        if [ "$emapper_types" -gt 0 ]; then
            echo "PASS: conduit_annotations.txt contains $emapper_types emapper-derived annotation type(s)"
        else
            echo "FAIL: conduit_annotations.txt contains no emapper-derived annotation types (protein ID mismatch?)"
            FAILED+=("uniprot_proteome_id_annotation_validation")
        fi
    fi
}
run_MAGs()                { run_integration_test MAGs; }
run_metaphlan() {
    # Seed the mock merged_profiles.txt into the run directory so downstream
    # rules are not pruned from the DAG when merge_profiles is omitted.
    local mock_src="$REPO_ROOT/experiments/integration_test/input/metaphlan/merged_profiles.txt"
    local mock_dst="$REPO_ROOT/experiments/integration_test/runs/metaphlan/metaphlan/merged_profiles.txt"
    mkdir -p "$(dirname "$mock_dst")"
    cp "$mock_src" "$mock_dst"

    # Skip MetaPhlAn database download and profiling steps.
    run_integration_test metaphlan \
        "--omit-from download_metaphlan_resources run_metaphlan merge_profiles"
}
run_peptidotyping() {
    if [ ! -f "$REPO_ROOT/resources/peptidotyping/sequences.tsv.lz4" ]; then
        echo "SKIP: peptidotyping — sequence index not found at resources/peptidotyping/"
        echo "      Run the build_sequence_index rule first to generate this resource."
        return
    fi
    run_integration_test peptidotyping
}

case "$METHOD" in
    ncbi_taxonomy_id)    run_ncbi_taxonomy_id ;;
    uniprot_proteome_id) run_uniprot_proteome_id ;;
    MAGs)                run_MAGs ;;
    metaphlan)           run_metaphlan ;;
    peptidotyping)       run_peptidotyping ;;
    all)
        run_ncbi_taxonomy_id
        run_uniprot_proteome_id
        run_MAGs
        run_metaphlan
        run_peptidotyping
        ;;
    *)
        echo "Unknown method: $METHOD"
        echo "Valid options: ncbi_taxonomy_id | uniprot_proteome_id | MAGs | metaphlan | peptidotyping | all"
        exit 1
        ;;
esac

# ── Results ───────────────────────────────────────────────────────────────────

echo "============================================================"
echo "Results: ${#PASSED[@]} passed, ${#FAILED[@]} failed"
if [ ${#PASSED[@]} -gt 0 ]; then
    echo "  PASSED: ${PASSED[*]}"
fi
if [ ${#FAILED[@]} -gt 0 ]; then
    echo "  FAILED: ${FAILED[*]}"
    exit 1
fi
