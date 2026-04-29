#!/bin/bash
# =============================================================================
# Create Test Sequences File
# =============================================================================
# Creates a smaller test file from the full sequences.tsv.lz4 file.
# Takes the first 10 million lines and saves it back as a compressed file.
# =============================================================================

set -euo pipefail

# Default values
INPUT_FILE="${1:-}"
OUTPUT_FILE="${2:-}"
SAMPLE_SIZE="${3:-10000000}"  # Default: 10 million lines

if [ -z "$INPUT_FILE" ] || [ -z "$OUTPUT_FILE" ]; then
    echo "Usage: $0 <input_sequences.tsv.lz4> <output_sequences.tsv.lz4> [sample_size]"
    echo ""
    echo "Arguments:"
    echo "  input_sequences.tsv.lz4  - Path to input sequences file (compressed)"
    echo "  output_sequences.tsv.lz4 - Path to output test file (compressed)"
    echo "  sample_size              - Number of lines to take (default: 10000000)"
    echo ""
    echo "Example:"
    echo "  $0 sequences.tsv.lz4 test_sequences.tsv.lz4 10000000"
    exit 1
fi

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file not found: $INPUT_FILE"
    exit 1
fi

echo "Creating test file from: $INPUT_FILE"
echo "Output file: $OUTPUT_FILE"
echo "Taking first $SAMPLE_SIZE lines..."
echo ""

# Decompress, take first N lines, recompress
lz4 -d -c "$INPUT_FILE" | head -n "$SAMPLE_SIZE" | lz4 - > "$OUTPUT_FILE"

# Get file sizes for comparison
INPUT_SIZE=$(du -h "$INPUT_FILE" | cut -f1)
OUTPUT_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)

echo ""
echo "Done!"
echo "Input size:  $INPUT_SIZE"
echo "Output size: $OUTPUT_SIZE"
echo ""
echo "To use this test file, update your config to point to: $OUTPUT_FILE"


