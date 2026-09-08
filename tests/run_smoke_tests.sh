#!/usr/bin/env bash
# Smoke runner: end-to-end pipeline tests against subset fixtures.
#
# The smoke layer sits between Tier 1 unit tests (no DBs, run in CI on every
# push) and Tier 3 full integration (real ~170 GB UMGAP + bakta + eggNOG DBs,
# manual). It points peptidotyping_resource_dir at tests/fixtures/peptidotyping_subset/
# — a ~17 MB curated subset covering the 10 ATCC pool species — so DIA-NN can
# actually run end-to-end in minutes instead of hours.
#
# Usage:
#   bash tests/run_smoke_tests.sh [method]
#
# Methods (one config each, no mode variants):
#   unipept_peptidotyping   peptidotyping smoke against the subset DB
#                            (4 expected pool families detected)
#
# Prerequisites:
#   - snakemake >= 8, apptainer/singularity
#   - Profile at profiles/<hostname>/config.yaml (auto-detected)
#   - For unipept_peptidotyping: the host profile must point taxonkit_db_dir at
#     a populated NCBI taxonkit DB (peptidotyping rules use it via container).
#
# Smoke runs vs the integration tests:
#   - Smoke uses subset peptidotyping fixture; integration uses the full ~170
#     GB UMGAP index.
#   - Smoke is local-only (not in CI; subset fixture is LFS-tracked but heavy).
#   - Both invoke the same Snakefile and rules; only the resource_dir overrides
#     differ. If a smoke run passes, the production code path is exercised.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
METHOD="${1:-all}"

PASSED=()
FAILED=()
SKIPPED=()

# ── Preflight ────────────────────────────────────────────────────────────────

if [ ! -f "$REPO_ROOT/tests/data/sample1.mzML" ] && [ ! -f "$REPO_ROOT/tests/data/sample1.raw" ]; then
    echo "ERROR: neither tests/data/sample1.mzML nor tests/data/sample1.raw is present."
    echo "Ensure Git LFS is installed and run: git lfs pull"
    exit 1
fi

# ── Host profile (carries machine-specific bind mounts + resource paths) ──────

PROFILE_ARGS=()
_host_profile="$REPO_ROOT/profiles/$(hostname)"
if [ -d "$_host_profile" ] && [ -f "$_host_profile/config.yaml" ]; then
    PROFILE_ARGS+=(--profile "$_host_profile")
fi

# ── DIA-NN availability gate ──────────────────────────────────────────────────
# Conduit no longer ships DIA-NN (licence: one backup copy, no sublicensing).
# A smoke run performs a real DIA-NN search, so skip cleanly rather than fail
# when no installation is reachable. See "Obtaining DIA-NN" in the README.
if [ -z "${CONDUIT_DIANN_PATH:-}" ] \
   && ! compgen -G "$REPO_ROOT/resources/diann/*/diann-linux" > /dev/null 2>&1 \
   && ! compgen -G "$REPO_ROOT/resources/diann/*.AppImage" > /dev/null 2>&1 \
   && ! grep -q 'diann_path' "$_host_profile/config.yaml" 2>/dev/null; then
    echo "SKIP: no DIA-NN installation found (set CONDUIT_DIANN_PATH, or unzip" >&2
    echo "      DIA-NN-<version>-Academia-Linux.zip into resources/diann/)." >&2
    echo "      https://github.com/vdemichev/DiaNN/releases" >&2
    exit 0
fi

# ── Profile-level resource keys that the smoke run still needs ────────────────
# The host profile passes these via `--config K=V ...`. When the smoke runner
# adds its own `--config peptidotyping_resource_dir=...`, snakemake REPLACES
# the entire profile --config block (does not merge), so we have to re-emit
# every key we want preserved. Read them from the profile and forward unchanged
# except peptidotyping_resource_dir, which we override per smoke fixture.
PROFILE_CONFIG_OVERRIDES=()
if [ ${#PROFILE_ARGS[@]} -gt 0 ]; then
    while IFS= read -r line; do
        # Lines look like "  - key=value" under the `config:` block; strip prefix.
        kv="${line#*- }"
        kv="${kv//[$'\t\r\n ']}"
        case "$kv" in
            peptidotyping_resource_dir=*) ;;   # smoke overrides this
            *=*) PROFILE_CONFIG_OVERRIDES+=("$kv") ;;
        esac
    done < <(grep -E "^\s*-\s+[a-zA-Z_]+=" "$_host_profile/config.yaml" 2>/dev/null || true)
fi

# ── Per-method runner ────────────────────────────────────────────────────────

run_smoke_test() {
    local method="$1"
    local target="$2"      # specific output file the smoke targets (NOT rule all)
    local config="$REPO_ROOT/tests/configs/smoke_${method}.yaml"

    if [ ! -f "$config" ]; then
        echo "SKIP: $method — no config at $config"
        SKIPPED+=("$method")
        return
    fi

    echo "============================================================"
    echo "Smoke test: $method  →  $target"
    echo "============================================================"

    local cores_args=()
    if [ ${#PROFILE_ARGS[@]} -eq 0 ]; then
        cores_args+=(--cores 4)
    fi

    local cmd=(
        snakemake
        --snakefile "$REPO_ROOT/Snakefile"
        --configfile "$config"
        --use-singularity
        "${PROFILE_ARGS[@]}"
        "${cores_args[@]}"
    )
    # Inject preserved profile keys + the smoke fixture override AS ONE
    # --config block (so all overrides come from a single, last-wins block).
    # Followed by `--` and the specific output file we want built — short of
    # the full rule all, which would cascade into the annotation pipeline.
    if [ ${#PROFILE_CONFIG_OVERRIDES[@]} -gt 0 ]; then
        cmd+=(--config)
        cmd+=("${PROFILE_CONFIG_OVERRIDES[@]}")
        cmd+=("peptidotyping_resource_dir=$REPO_ROOT/tests/fixtures/peptidotyping_subset/")
    fi
    cmd+=(-- "$target")

    if "${cmd[@]}" 2>&1; then
        echo "PASS: $method"
        PASSED+=("$method")
    else
        echo "FAIL: $method"
        FAILED+=("$method")
    fi
    echo
}

# ── Dispatch ─────────────────────────────────────────────────────────────────

run_unipept_peptidotyping() {
    # Target the family-detection output rather than rule all. Going further
    # exercises get_uniprot_proteome_ids, which crashes on empty species/strain
    # input via conduitR::get_better_proteome_ids — a separate upstream bug,
    # not what smoke is meant to validate. The family-detection output is
    # what proves the peptidotyping decision logic works end-to-end.
    local target="experiments/integration_test/runs/smoke_unipept_peptidotyping/database_resources/peptidotyping/detected_family_taxa_ids.txt"
    run_smoke_test unipept_peptidotyping "$target"
}

case "$METHOD" in
    unipept_peptidotyping) run_unipept_peptidotyping ;;
    all)
        run_unipept_peptidotyping
        ;;
    *)
        echo "ERROR: unknown method '$METHOD'"
        echo "Valid options: unipept_peptidotyping | all"
        exit 1
        ;;
esac

# ── Summary ──────────────────────────────────────────────────────────────────

echo "============================================================"
echo "Results: ${#PASSED[@]} passed, ${#FAILED[@]} failed, ${#SKIPPED[@]} skipped"
if [ ${#PASSED[@]} -gt 0 ]; then
    echo "  PASSED: ${PASSED[*]}"
fi
if [ ${#SKIPPED[@]} -gt 0 ]; then
    echo "  SKIPPED: ${SKIPPED[*]}"
fi
if [ ${#FAILED[@]} -gt 0 ]; then
    echo "  FAILED: ${FAILED[*]}"
    exit 1
fi
