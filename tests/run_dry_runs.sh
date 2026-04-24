#!/usr/bin/env bash
# Run Snakemake dry-runs for all search_space_method configurations.
# Activate the snakemake conda environment before running:
#   conda activate snakemake && bash tests/run_dry_runs.sh
#
# Env-var overrides (see tests/run_integration_tests.sh for full docs):
#   PEPTIDOTYPING_RESOURCE_DIR, TAXONKIT_DB_DIR — off-repo paths for heavy resources
#   SINGULARITY_BIND_PATHS — apptainer bind paths (not strictly needed for dry-runs,
#     but accepted for symmetry with the integration runner)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIGS_DIR="$REPO_ROOT/tests/configs"

PASSED=()
FAILED=()

# Collect --config overrides from env vars (same contract as run_integration_tests.sh).
_config_overrides=()
if [ -n "${PEPTIDOTYPING_RESOURCE_DIR:-}" ]; then
    _config_overrides+=("peptidotyping_resource_dir=$PEPTIDOTYPING_RESOURCE_DIR")
fi
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

# ── macOS AppleDouble dropout: .yaml files starting with "._" are junk metadata, skip them.

run_dry_run() {
    local config_file="$1"
    local method
    method="$(basename "$config_file" .yaml)"

    echo "------------------------------------------------------------"
    echo "Dry-run: $method"
    echo "------------------------------------------------------------"

    if snakemake \
        --snakefile "$REPO_ROOT/Snakefile" \
        --configfile "$config_file" \
        --dry-run \
        --cores 1 \
        --quiet \
        "${EXTRA_CONFIG_ARGS[@]}" \
        2>&1; then
        echo "PASS: $method"
        PASSED+=("$method")
    else
        echo "FAIL: $method"
        FAILED+=("$method")
    fi
    echo
}

for config in "$CONFIGS_DIR"/*.yaml; do
    # Skip macOS resource-fork files (._*).
    [[ "$(basename "$config")" == ._* ]] && continue
    run_dry_run "$config"
done

echo "============================================================"
echo "Results: ${#PASSED[@]} passed, ${#FAILED[@]} failed"
if [ ${#PASSED[@]} -gt 0 ]; then
    echo "  PASSED: ${PASSED[*]}"
fi
if [ ${#FAILED[@]} -gt 0 ]; then
    echo "  FAILED: ${FAILED[*]}"
    exit 1
fi
