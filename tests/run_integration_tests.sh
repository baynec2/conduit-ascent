#!/usr/bin/env bash
# Run end-to-end integration tests for one or all search_space_method configurations.
#
# Usage:
#   bash tests/run_integration_tests.sh [method]
#
#   method: preflight | ncbi_taxonomy_id | uniprot_proteome_id | MAGs | metaphlan
#           | unipept_peptidotyping | all
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
# Per-machine configuration lives in profiles/$(hostname)/config.yaml.
# The runner auto-detects a directory matching the current hostname and
# passes --profile to snakemake. Put this machine's bind mounts, resource
# path overrides, and cores there. See profiles/nanopore-catalyst/ for
# the reference layout; copy that directory name → your hostname and
# edit paths to match your storage.
#
# Without a host profile, the portable base-config defaults apply — fine
# for ncbi_taxonomy_id / uniprot_proteome_id / MAGs / metaphlan / hapid,
# insufficient for unipept_peptidotyping family (those methods require the ~170 GB
# UMGAP index, which lives off-repo on most machines).
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
# The unipept_peptidotyping method requires a pre-built sequence index at:
#   resources/peptidotyping/  (or wherever the host profile's
#   peptidotyping_resource_dir points)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
METHOD="${1:-all}"

PASSED=()
FAILED=()

# ── Preflight checks ──────────────────────────────────────────────────────────

# Require at least one test spectra file — either Thermo .raw or open-format
# .mzML. Both are LFS-tracked; the current default is sample1.mzML.
if [ ! -f "$REPO_ROOT/tests/data/sample1.mzML" ] && [ ! -f "$REPO_ROOT/tests/data/sample1.raw" ]; then
    echo "ERROR: neither tests/data/sample1.mzML nor tests/data/sample1.raw is present."
    echo "Ensure Git LFS is installed and run: git lfs pull"
    exit 1
fi

# ── Host-matched Snakemake profile ────────────────────────────────────────────

# If profiles/$(hostname)/config.yaml exists, apply it. The profile carries
# machine-specific settings (bind mounts, resource overrides, cores).
PROFILE_ARGS=()
_host_profile="$REPO_ROOT/profiles/$(hostname)"
if [ -d "$_host_profile" ] && [ -f "$_host_profile/config.yaml" ]; then
    PROFILE_ARGS+=(--profile "$_host_profile")
fi

# ── Peptidotyping resource pre-check ──────────────────────────────────────────
# Skip unipept_peptidotyping-family tests with a helpful message when the index is
# unavailable. With a host profile active, trust the profile (it's responsible
# for pointing at the right path); otherwise check the portable default.
if [ ${#PROFILE_ARGS[@]} -eq 0 ]; then
    PEPTIDOTYPING_INDEX="$REPO_ROOT/resources/peptidotyping/sequences.tsv.lz4"
else
    PEPTIDOTYPING_INDEX=""   # bypass pre-check; let snakemake resolve
fi

# ── Runner ────────────────────────────────────────────────────────────────────

run_integration_test() {
    local method="$1"
    local extra_flags_str="${2:-}"
    local config="$REPO_ROOT/tests/configs/integration_test_${method}.yaml"

    echo "============================================================"
    echo "Integration test: $method"
    echo "============================================================"

    # With a host profile active, let it set cores (typically 20 for a
    # workstation). Without a profile (e.g., CI on a tiny runner), fall back
    # to 4 cores so snakemake has something to work with.
    local cores_args=()
    if [ ${#PROFILE_ARGS[@]} -eq 0 ]; then
        cores_args+=(--cores 4)
    fi

    local snakemake_cmd=(
        snakemake
        --snakefile "$REPO_ROOT/Snakefile"
        --configfile "$config"
        --use-singularity
        "${PROFILE_ARGS[@]}"
        "${cores_args[@]}"
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
    # The integration_test sample1.fastq.gz is a 1k-read subset of an E. coli
    # sequencing run (sourced from experiments/s76/input/fastq_files/76.fastq.gz)
    # — small enough for fast iteration, large enough for MetaPhlAn to hit
    # E. coli marker genes reliably. The MetaPhlAn DB is expected at
    # config["metaphlan_database_dir"]; on a host with the DB already present
    # the install rule is a no-op (only the sentinel mpa_latest file is
    # checked). No --omit-from games — the full DAG runs.
    run_integration_test metaphlan
}
run_unipept_peptidotyping() {
    if [ -n "$PEPTIDOTYPING_INDEX" ] && [ ! -f "$PEPTIDOTYPING_INDEX" ]; then
        echo "SKIP: unipept_peptidotyping — sequence index not found at resources/peptidotyping/"
        echo "      Run the build_sequence_index rule first to generate this resource."
        return
    fi
    run_integration_test unipept_peptidotyping
}
run_genome_peptidotyping() {
    # Shares the unipept_peptidotyping sequence-index prerequisite.
    if [ -n "$PEPTIDOTYPING_INDEX" ] && [ ! -f "$PEPTIDOTYPING_INDEX" ]; then
        echo "SKIP: genome_peptidotyping — unipept_peptidotyping sequence index not found."
        return
    fi
    run_integration_test genome_peptidotyping
}
run_unipept_hapid() {
    # Shares the unipept_peptidotyping sequence-index prerequisite.
    if [ -n "$PEPTIDOTYPING_INDEX" ] && [ ! -f "$PEPTIDOTYPING_INDEX" ]; then
        echo "SKIP: unipept_hapid — unipept_peptidotyping sequence index not found."
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
    unipept_peptidotyping) run_unipept_peptidotyping ;;
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
        run_unipept_peptidotyping
        run_genome_peptidotyping
        run_unipept_hapid
        run_hapid
        run_mgnify_MAGs
        run_mgnify_hapid
        ;;
    *)
        echo "Unknown method: $METHOD"
        echo "Valid options: preflight | ncbi_taxonomy_id | uniprot_proteome_id | MAGs | metaphlan"
        echo "             | unipept_peptidotyping | genome_peptidotyping | unipept_hapid | hapid"
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
