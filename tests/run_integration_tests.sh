#!/usr/bin/env bash
# Run end-to-end integration tests for one or all search_space_method configurations.
#
# Usage:
#   bash tests/run_integration_tests.sh [method]
#
#   method: preflight | ncbi_taxonomy_id | uniprot_proteome_id | MAGs | metaphlan
#           | peptidotyping | all
#   Defaults to "all" if not specified.
#
# The "preflight" target runs DIA-NN against tests/data/sample1.raw with a trivial
# FASTA to confirm the raw file is readable. When invoked via "all", preflight runs
# first and aborts the suite on failure (no point running downstream tests if the
# raw file cannot be read).
#
# Prerequisites:
#   - snakemake >= 8
#   - apptainer/singularity
#   - tests/data/sample1.raw (tracked via Git LFS)
#
# Env-var overrides (optional, per-machine; test configs stay portable):
#   PEPTIDOTYPING_RESOURCE_DIR   if the UMGAP index lives off-repo (e.g., HDD),
#                                set this and the runner will pass
#                                --config peptidotyping_resource_dir=$VALUE
#                                to snakemake. Takes precedence over the config-file value.
#   TAXONKIT_DB_DIR              location of the taxonkit NCBI taxonomy DB (required by
#                                peptidotyping + genome_peptidotyping). If
#                                PEPTIDOTYPING_RESOURCE_DIR is set and this is not,
#                                defaults to $PEPTIDOTYPING_RESOURCE_DIR/taxonkit per
#                                the CLAUDE.md convention.
#   SINGULARITY_BIND_PATHS       comma- or colon-separated list of paths to bind
#                                into apptainer containers (e.g., /home/nanopore-catalyst/HDD).
#                                Required whenever an overridden resource dir is outside
#                                the repo tree, since apptainer does not see external
#                                paths by default.
#
# The metaphlan method mocks the MetaPhlAn database + profiling steps via --omit-from.
# A pre-committed mock merged_profiles.txt is seeded into the run directory so that
# downstream rules (call_ncbi_taxa_ids and beyond) are not also pruned from the DAG.
#
# The MAGs method requires a FASTA file at:
#   experiments/integration_test/input/MAG_files/atcc_25922.fasta
#
# All methods read MS spectra from experiments/<exp>/input/ms_files/ — .raw and
# .mzML are both accepted by DIA-NN. For integration tests this directory
# contains a symlink to tests/data/sample1.mzML.
#
# The peptidotyping method requires a pre-built sequence index at:
#   resources/peptidotyping/  (or wherever PEPTIDOTYPING_RESOURCE_DIR points)

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

# ── Env-var overrides → snakemake args ────────────────────────────────────────

# Collect all config overrides into a single --config invocation, since snakemake's
# --config uses `store` semantics (a second --config overwrites the first).
_config_overrides=()
if [ -n "${PEPTIDOTYPING_RESOURCE_DIR:-}" ]; then
    _config_overrides+=("peptidotyping_resource_dir=$PEPTIDOTYPING_RESOURCE_DIR")
fi
# Taxonkit DB: explicit env var wins; else derive from PEPTIDOTYPING_RESOURCE_DIR/taxonkit
# per the convention documented in CLAUDE.md.
_taxonkit_db="${TAXONKIT_DB_DIR:-}"
if [ -z "$_taxonkit_db" ] && [ -n "${PEPTIDOTYPING_RESOURCE_DIR:-}" ]; then
    _taxonkit_db="${PEPTIDOTYPING_RESOURCE_DIR%/}/taxonkit"
fi
if [ -n "$_taxonkit_db" ]; then
    _config_overrides+=("taxonkit_db_dir=$_taxonkit_db")
fi

EXTRA_CONFIG_ARGS=()
if [ ${#_config_overrides[@]} -gt 0 ]; then
    EXTRA_CONFIG_ARGS+=(--config "${_config_overrides[@]}")
fi

SINGULARITY_ARGS=()
if [ -n "${SINGULARITY_BIND_PATHS:-}" ]; then
    # Normalize separators (accept comma or colon) → space.
    bind_arg="--bind ${SINGULARITY_BIND_PATHS//[,:]/,}"
    SINGULARITY_ARGS+=(--singularity-args "$bind_arg")
fi

# ── Peptidotyping resource resolution ─────────────────────────────────────────
# Used by method dispatch functions to skip when the index is unavailable.
_pep_root="${PEPTIDOTYPING_RESOURCE_DIR:-$REPO_ROOT/resources/peptidotyping}"
PEPTIDOTYPING_INDEX="${_pep_root%/}/sequences.tsv.lz4"

# ── Runner ────────────────────────────────────────────────────────────────────

run_integration_test() {
    local method="$1"
    local extra_flags_str="${2:-}"
    local config="$REPO_ROOT/tests/configs/integration_test_${method}.yaml"

    echo "============================================================"
    echo "Integration test: $method"
    echo "============================================================"

    local snakemake_cmd=(
        snakemake
        --snakefile "$REPO_ROOT/Snakefile"
        --configfile "$config"
        --use-singularity
        --cores 4
        "${SINGULARITY_ARGS[@]}"
        "${EXTRA_CONFIG_ARGS[@]}"
    )
    # extra_flags_str is a string (e.g., "--omit-from foo bar"). Word-split it.
    if [ -n "$extra_flags_str" ]; then
        # shellcheck disable=SC2206
        local extra_flags=($extra_flags_str)
        snakemake_cmd+=("${extra_flags[@]}")
    fi

    if "${snakemake_cmd[@]}" 2>&1; then
        echo "PASS: $method"
        PASSED+=("$method")
    else
        echo "FAIL: $method"
        FAILED+=("$method")
    fi
    echo
}

# ── Method dispatch ───────────────────────────────────────────────────────────

run_preflight() {
    # Preflight exercises DIA-NN's library-search code path (--dir + --lib) —
    # the same path run_diann uses in production. An earlier version only
    # exercised --fasta-search + --f, which passed even on a file that broke
    # the real workflow (see memory/feedback_preflight_must_match_target_codepath.md).
    #
    # Steps:
    #   1. Discover the test spectra file (sample1.mzML preferred; sample1.raw fallback)
    #   2. Predict a tiny library from a one-protein FASTA
    #   3. Run DIA-NN library search with --dir pointing at the test data
    #   4. Parse the scan-count line to confirm the file loaded
    echo "============================================================"
    echo "Preflight: DIA-NN --dir + --lib read test on tests/data/"
    echo "============================================================"
    local diann_image="docker://baynec2/diann2.1.0:alpha"

    # Pick an input — prefer mzML (open format, portable); fall back to .raw.
    local spectra_file
    if [ -f "$REPO_ROOT/tests/data/sample1.mzML" ]; then
        spectra_file="sample1.mzML"
    elif [ -f "$REPO_ROOT/tests/data/sample1.raw" ]; then
        spectra_file="sample1.raw"
    else
        echo "FAIL: preflight — neither tests/data/sample1.mzML nor sample1.raw found"
        FAILED+=("preflight")
        return
    fi

    local tmpdir
    tmpdir="$(mktemp -d)"
    # Isolate the single test file into its own dir — DIA-NN --dir would
    # otherwise pick up any other files (sidecars, stale outputs) alongside.
    mkdir -p "$tmpdir/raw"
    ln -sf "$REPO_ROOT/tests/data/$spectra_file" "$tmpdir/raw/$spectra_file"
    printf '>preflight_test_protein\nMSSSTPPAQKGRGTQGRGRVSLGEGHLGFTTKPSIGGELHTEEA\n' \
        > "$tmpdir/tiny.fasta"

    local predict_log="$tmpdir/predict.log"
    local search_log="$tmpdir/search.log"

    # Step 1: predict a tiny library from the tiny FASTA.
    apptainer exec --bind "$tmpdir:/work" "$diann_image" \
        diann \
            --fasta /work/tiny.fasta \
            --fasta-search --gen-spec-lib --predictor \
            --cut "K*,R*" --missed-cleavages 0 \
            --unimod4 --mass-acc 10 --mass-acc-ms1 4 \
            --min-pep-len 7 --max-pep-len 30 \
            --out-lib /work/tiny \
            --threads 4 \
            > "$predict_log" 2>&1 || true
    if [ ! -f "$tmpdir/tiny.predicted.speclib" ]; then
        echo "FAIL: preflight — library prediction step did not produce /tmp/.../tiny.predicted.speclib"
        echo "Last 30 lines of predict log:"
        tail -30 "$predict_log"
        FAILED+=("preflight")
        return
    fi

    # Step 2: library-based search against /work/raw (the real path run_diann uses).
    apptainer exec --bind "$tmpdir:/work" "$diann_image" \
        diann \
            --dir /work/raw \
            --lib /work/tiny.predicted.speclib \
            --fasta /work/tiny.fasta \
            --qvalue 1 \
            --mass-acc 10 --mass-acc-ms1 4 \
            --out /work/search_report \
            --threads 4 \
            > "$search_log" 2>&1 || true

    local scan_line
    scan_line=$(grep -E "[0-9]+ MS1 and [0-9]+ MS2 scans" "$search_log" | head -1) || true
    if [ -z "$scan_line" ]; then
        echo "FAIL: preflight — DIA-NN library search did not report scan counts for $spectra_file"
        echo "Last 40 lines of search log ($search_log):"
        tail -40 "$search_log"
        FAILED+=("preflight")
        return
    fi
    # Extract MS1/MS2 counts by field, not regex — "MS1"/"MS2" contain digits that
    # a naive [0-9]+ match would pick up alongside the intended scan counts.
    local ms1 ms2
    ms1=$(echo "$scan_line" | awk '{for(i=1;i<NF;i++) if($(i+1)=="MS1") {print $i; exit}}')
    ms2=$(echo "$scan_line" | awk '{for(i=1;i<NF;i++) if($(i+1)=="MS2") {print $i; exit}}')
    if [ "${ms1:-0}" -gt 0 ] && [ "${ms2:-0}" -gt 0 ]; then
        echo "PASS: preflight — $ms1 MS1 and $ms2 MS2 scans read from $spectra_file (--dir + --lib mode)"
        PASSED+=("preflight")
        rm -rf "$tmpdir"
    else
        echo "FAIL: preflight — zero scans read ($ms1 MS1, $ms2 MS2)"
        echo "Log preserved at: $search_log"
        FAILED+=("preflight")
    fi
    echo
}

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
    if [ ! -f "$PEPTIDOTYPING_INDEX" ]; then
        echo "SKIP: peptidotyping — sequence index not found at resources/peptidotyping/"
        echo "      Run the build_sequence_index rule first to generate this resource."
        return
    fi
    run_integration_test peptidotyping
}
run_genome_peptidotyping() {
    # Shares the peptidotyping sequence-index prerequisite.
    if [ ! -f "$PEPTIDOTYPING_INDEX" ]; then
        echo "SKIP: genome_peptidotyping — peptidotyping sequence index not found."
        return
    fi
    run_integration_test genome_peptidotyping
}
run_unipept_hapid() {
    # Shares the peptidotyping sequence-index prerequisite.
    if [ ! -f "$PEPTIDOTYPING_INDEX" ]; then
        echo "SKIP: unipept_hapid — peptidotyping sequence index not found."
        return
    fi
    run_integration_test unipept_hapid
}
run_hapid() {
    if [ ! -f "$REPO_ROOT/resources/hapid/ribP_elonF_profiles_refined_manually.hmm" ]; then
        echo "SKIP: hapid — HAPiID HMM profiles not found at resources/hapid/."
        return
    fi
    run_integration_test hapid
}
run_mgnify_MAGs() {
    # Exercises the MGnify download path; mgnify_max_genomes: 3 keeps it cheap.
    run_integration_test mgnify_MAGs
}
run_mgnify_hapid() {
    if [ ! -f "$REPO_ROOT/resources/hapid/ribP_elonF_profiles_refined_manually.hmm" ]; then
        echo "SKIP: mgnify_hapid — HAPiID HMM profiles not found at resources/hapid/."
        return
    fi
    run_integration_test mgnify_hapid
}

case "$METHOD" in
    preflight)           run_preflight ;;
    ncbi_taxonomy_id)    run_ncbi_taxonomy_id ;;
    uniprot_proteome_id) run_uniprot_proteome_id ;;
    MAGs)                run_MAGs ;;
    metaphlan)           run_metaphlan ;;
    peptidotyping)       run_peptidotyping ;;
    genome_peptidotyping) run_genome_peptidotyping ;;
    unipept_hapid)       run_unipept_hapid ;;
    hapid)               run_hapid ;;
    mgnify_MAGs)         run_mgnify_MAGs ;;
    mgnify_hapid)        run_mgnify_hapid ;;
    all)
        run_preflight
        if [ ${#FAILED[@]} -gt 0 ]; then
            echo "ERROR: preflight failed — aborting integration tests."
            echo "  FAILED: ${FAILED[*]}"
            exit 1
        fi
        run_ncbi_taxonomy_id
        run_uniprot_proteome_id
        run_MAGs
        run_metaphlan
        run_peptidotyping
        run_genome_peptidotyping
        run_unipept_hapid
        run_hapid
        run_mgnify_MAGs
        run_mgnify_hapid
        ;;
    *)
        echo "Unknown method: $METHOD"
        echo "Valid options: preflight | ncbi_taxonomy_id | uniprot_proteome_id | MAGs | metaphlan"
        echo "             | peptidotyping | genome_peptidotyping | unipept_hapid | hapid"
        echo "             | mgnify_MAGs | mgnify_hapid | all"
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
