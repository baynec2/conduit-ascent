#!/usr/bin/env bash
# Generate a small mzML test file from a Thermo .raw for integration tests.
#
# Why not ship a .raw directly? Some .raw files lose instrument-index
# metadata during filtering / subsetting — DIA-NN's library-search code path
# needs that metadata and fails with "ThermoRaw exception: Instrument index
# not available". mzML is an open format with well-specified headers, so
# converting once and committing the mzML avoids that whole class of
# problems. See memory/feedback_preflight_must_match_target_codepath.md
# for the incident context.
#
# Usage:
#   bash tests/scripts/generate_test_mzml.sh <input.raw> [<scan_start> <scan_end>]
#
# Example:
#   bash tests/scripts/generate_test_mzml.sh \
#       /path/to/known_good_calibration.raw 1000 2000
#
# Requirements:
#   - apptainer/singularity (for ThermoRawFileParser container)
#   - msconvert on PATH  (for mzML → mzML scan-range subsetting)
#
# Output:
#   tests/data/sample1.mzML   (subset if scan range given, else full)
#   tests/data/sample1_full.mzML.size.txt  (size record, diagnostic)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT_DIR="$REPO_ROOT/tests/data"

INPUT_RAW="${1:-}"
SCAN_START="${2:-}"
SCAN_END="${3:-}"

if [ -z "$INPUT_RAW" ] || [ ! -f "$INPUT_RAW" ]; then
    echo "ERROR: pass a valid Thermo .raw file as the first argument" >&2
    echo "Usage: $0 <input.raw> [<scan_start> <scan_end>]" >&2
    exit 1
fi

if ! command -v apptainer >/dev/null && ! command -v singularity >/dev/null; then
    echo "ERROR: apptainer or singularity is required for ThermoRawFileParser" >&2
    exit 1
fi

if ! command -v msconvert >/dev/null; then
    echo "ERROR: msconvert not on PATH; install ProteoWizard (mzML → mzML subsetting)" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

# ── Step 1: Thermo .raw → mzML via ThermoRawFileParser ────────────────────────

FULL_MZML="$OUT_DIR/sample1_full.mzML"
echo "→ Converting $INPUT_RAW → $FULL_MZML via ThermoRawFileParser"

INPUT_DIR="$(cd "$(dirname "$INPUT_RAW")" && pwd)"
INPUT_BASE="$(basename "$INPUT_RAW")"

# --format=2 emits INDEXED mzML. The index at the end of the file is not
# optional for DIA-NN: without it DIA-NN reports "No MS2 spectra: aborting"
# even when the file has valid MS2 scans. Non-indexed mzML (--format=1) is
# smaller but unusable here. This is the subtle bug that caught us on the
# first pass of Phase 4 smoke testing.
apptainer exec \
    --bind "$INPUT_DIR:/in" \
    --bind "$OUT_DIR:/out" \
    docker://quay.io/biocontainers/thermorawfileparser:1.4.3--ha8f3691_0 \
    ThermoRawFileParser.sh \
        --input "/in/$INPUT_BASE" \
        --output_file "/out/sample1_full.mzML" \
        --format=2

echo "  $(du -h "$FULL_MZML" | awk '{print $1}')  full mzML"

# ── Step 2: Optional scan-range subset via msconvert ──────────────────────────

FINAL_MZML="$OUT_DIR/sample1.mzML"

if [ -n "$SCAN_START" ] && [ -n "$SCAN_END" ]; then
    echo "→ Subsetting scans $SCAN_START..$SCAN_END via msconvert → $FINAL_MZML"
    # msconvert writes into --outdir with --outfile; it refuses to overwrite
    # unless the target isn't the same name as the input, so we use a
    # temp name then rename.
    msconvert "$FULL_MZML" \
        --filter "scanNumber [$SCAN_START,$SCAN_END]" \
        --mzML --zlib \
        --outdir "$OUT_DIR" \
        --outfile "sample1.mzML"
else
    echo "→ No scan range given; using full mzML as tests/data/sample1.mzML"
    cp "$FULL_MZML" "$FINAL_MZML"
fi

echo "  $(du -h "$FINAL_MZML" | awk '{print $1}')  final mzML"

# ── Sanity diagnostics ────────────────────────────────────────────────────────

# grep -c prints the count to stdout then exits 1 on zero matches; disable
# `set -e` briefly so the 0-match case doesn't abort the script. Note that
# mzML ordering puts value= before name= on the ms-level cvParam.
set +e
MS1_COUNT=$(grep -c 'value="1" name="ms level"' "$FINAL_MZML")
MS2_COUNT=$(grep -c 'value="2" name="ms level"' "$FINAL_MZML")
set -e
echo "  MS1 spectra: $MS1_COUNT"
echo "  MS2 spectra: $MS2_COUNT"

if [ "$MS1_COUNT" -eq 0 ] && [ "$MS2_COUNT" -eq 0 ]; then
    echo "WARNING: no MS1/MS2 spectra visible in grep — verify the mzML is valid" >&2
fi

echo ""
echo "Done. Next:"
echo "  - Place a symlink at experiments/<exp>/input/ms_files/sample1.mzML"
echo "    pointing to tests/data/sample1.mzML."
echo "  - Confirm the Snakefile's RAW_FILEPATHS glob picks up *.mzML too."
echo "  - Update tests/run_integration_tests.sh preflight to use --dir + --lib mode"
echo "    so we catch library-search path regressions."
