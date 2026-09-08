#!/usr/bin/env bash
# Run Snakemake dry-runs for all search_space_method configurations.
# Activate the snakemake conda environment before running:
#   conda activate snakemake && bash tests/run_dry_runs.sh
#
# Per-machine config via profiles/$(hostname)/config.yaml (auto-picked up).
# See tests/run_integration_tests.sh for the full contract.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIGS_DIR="$REPO_ROOT/tests/configs"

PASSED=()
FAILED=()

# Host-matched profile, if present.
PROFILE_ARGS=()
_host_profile="$REPO_ROOT/profiles/$(hostname)"
if [ -d "$_host_profile" ] && [ -f "$_host_profile/config.yaml" ]; then
    PROFILE_ARGS+=(--profile "$_host_profile")
fi

# No DIA-NN gate here, unlike run_integration_tests.sh / run_smoke_tests.sh: a
# dry run never executes DIA-NN, and the Snakefile's preflight downgrades a
# missing installation to a warning under --dry-run so the DAG still builds.
# That is what lets these run in CI, which cannot hold the binary.

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
        "${PROFILE_ARGS[@]}" \
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
