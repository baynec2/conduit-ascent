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
#
# The build is split into two rules because its halves fail for unrelated
# reasons and each takes hours:
#
#   download_uniprot_peptides  network-bound; downloads SwissProt+TrEMBL and
#                              digests them into peptides-out.tsv.lz4 (~9.5 h)
#   build_sequence_index       disk-bound; sorts that table and derives the
#                              LCA/FA tables and sequences.tsv.lz4
#
# As one rule, an FTP hiccup or a failed sort discarded the whole thing. Split,
# Snakemake resumes at the boundary and `retries` can target the flaky half.
################################################################################

import os

# Where the build's scratch lives. Defaults under the resource dir; override it
# when that path is on a slow filesystem. On WSL the Windows drives are exposed
# over 9p with msize=65536 (64 KB per round trip), and GNU sort writes its spill
# files in 4 KB units -- filling 1/16 of every message. Measured: ~11 MB/s on
# /mnt/c versus ~4 GB/s on the VM's own ext4 disk, which is the difference
# between hours and days. Point this at a local disk on such hosts.
PT_RES = config["peptidotyping_resource_dir"]
PT_TMP = config.get("peptidotyping_temp_dir", os.path.join(PT_RES, "temp"))

UMGAP_SH = ("modules/search_space/unipept_peptidotyping/scripts/"
            "unipept-database/scripts/generate_umgap_tables.sh")

# Shared prologue. CARGO_HOME is redirected because the Dockerfile pins it to a
# path inside the image, which Apptainer mounts read-only; the runtime
# `cargo build --release` of the rust-utils helpers needs somewhere writable.
# A stable location under the resource dir (not temp, which gets cleaned) keeps
# the crate downloads across runs. RUSTUP_HOME is untouched, so the rustup cargo
# proxy still finds its toolchain.
_PROLOGUE = r"""
        set -euo pipefail
        mkdir -p {params.outdir} {params.temp_outdir} $(dirname {log})
        export TMPDIR=$(realpath {params.temp_outdir})
        export CARGO_HOME=$(realpath {params.outdir})/cargo
        mkdir -p "$CARGO_HOME"
"""


# This is needed to generate the file containing all peptides in TREMBL and
# SWISSPROT and their LCAs. See https://github.com/unipept/unipept-database/issues/75
rule download_uniprot_peptides:
    output:
        # temp(): removed once build_sequence_index succeeds, but KEPT if it
        # fails -- so a failed index build resumes from here instead of
        # re-downloading SwissProt+TrEMBL. Lives directly in the temp dir, not
        # in its unipept_temp/ subdirectory, which the script's EXIT trap wipes.
        peptides = temp(os.path.join(PT_TMP, "peptides-out.tsv.lz4")),
        taxons   = os.path.join(PT_RES, "taxons.tsv.lz4"),
        lineages = os.path.join(PT_RES, "lineages.tsv.lz4"),
        entries  = os.path.join(PT_RES, "uniprot_entries.tsv.lz4"),
        relnotes = os.path.join(PT_RES, "relnotes.txt")
    params:
        outdir      = PT_RES,
        temp_outdir = PT_TMP,
        script      = UMGAP_SH
    benchmark:
        os.path.join(PT_RES, "benchmarks/download_uniprot_peptides.tsv")
    log:
        os.path.join(PT_RES, "logs/download_uniprot_peptides.log")
    container:
        config["containers"]["umgap"]
    shell:
        _PROLOGUE + r"""
        # Record the UniProt release alongside the data it describes.
        curl -L -o {output.relnotes} \
          https://ftp.uniprot.org/pub/databases/uniprot/relnotes.txt

        ONLY=download {params.script} tryptic \
          --output-dir {params.outdir} \
          --database-sources swissprot,trembl \
          --temp-dir {params.temp_outdir} \
          --min-peptide-length 5 \
          --max-peptide-length 50 \
          >> {log} 2>&1
        """


rule build_sequence_index:
    input:
        peptides = rules.download_uniprot_peptides.output.peptides,
        taxons   = rules.download_uniprot_peptides.output.taxons
    output:
        sequences = os.path.join(PT_RES, "sequences.tsv.lz4")
    params:
        outdir      = PT_RES,
        temp_outdir = PT_TMP,
        script      = UMGAP_SH
    benchmark:
        os.path.join(PT_RES, "benchmarks/build_sequence_index.tsv")
    log:
        os.path.join(PT_RES, "logs/build_sequence_index.log")
    container:
        config["containers"]["umgap"]
    shell:
        _PROLOGUE + r"""
        ONLY=index {params.script} tryptic \
          --output-dir {params.outdir} \
          --temp-dir {params.temp_outdir} \
          >> {log} 2>&1
        """
