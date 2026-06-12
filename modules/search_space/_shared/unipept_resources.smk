################################################################################
# Shared Unipept resources
################################################################################
# The UMGAP-derived sequence index (sequences.tsv.lz4 + taxons.tsv.lz4) is
# consumed by both the peptidotyping and unipept_hapid search_space methods.
# This module is the canonical location for the build rule so the two methods
# can coexist in the same workflow DAG without rule-name collisions.
#
# Include via `use rule * from shared_unipept_resources` in the main Snakefile
# wherever either peptidotyping or unipept_hapid is active.
################################################################################

import os

# This is needed to generate the file containing all peptides in TREMBL and
# SWISSPROT and their LCAs. See https://github.com/unipept/unipept-database/issues/75
rule build_sequence_index:
    output:
        sequences = os.path.join(config["peptidotyping_resource_dir"],"sequences.tsv.lz4"),
        taxons    = os.path.join(config["peptidotyping_resource_dir"],"taxons.tsv.lz4")
    params:
        outdir      = config["peptidotyping_resource_dir"],
        temp_outdir = os.path.join(config["peptidotyping_resource_dir"],"temp")
    benchmark:
        os.path.join(config["peptidotyping_resource_dir"],"benchmarks/build_sequence_index.tsv")
    log:
        os.path.join(config["peptidotyping_resource_dir"],"logs/build_sequence_index.log")
    container:
        config["containers"]["umgap"]
    shell:
        r"""
        set -euo pipefail

        # Ensure directories exist
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})
        mkdir -p {params.temp_outdir}

        # Download UniProt release notes
        curl -L \
          -o {params.outdir}/relnotes.txt \
          https://ftp.uniprot.org/pub/databases/uniprot/relnotes.txt

        # Setting temp dir (must be absolute so cargo resolves it correctly)
        export TMPDIR=$(realpath {params.temp_outdir})

        # Build UMGAP peptidotyping tables
        modules/search_space/unipept_peptidotyping/scripts/unipept-database/scripts/generate_umgap_tables.sh tryptic \
          --output-dir {params.outdir} \
          --database-sources swissprot,trembl \
          --temp-dir {params.temp_outdir} \
          --min-peptide-length 5 \
          --max-peptide-length 50 \
          >> {log} 2>&1
        """
