import os

# ==============================================================================
# Paths
# ==============================================================================
EXPERIMENT_DIR = config["experiment_dir"]
RUN_DIR        = config["run_dir"]
MAG_DIR        = os.path.join(EXPERIMENT_DIR, "input/MAG_files")
MGNIFY_OUT     = os.path.join(RUN_DIR, "genome_download/mgnify")

MGNIFY_FTP_BASE = config.get("mgnify_ftp_base",
    "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes")
MGNIFY_CATALOG  = config.get("mgnify_catalog", "")

# Globalized cache for content-identical files: genomes and the raw
# metadata.tsv live under a per-catalog directory keyed by the catalog
# string (e.g., "human-gut/v2.0.2" → "human-gut_v2.0.2") so multiple
# experiments / runs using the same catalog share one copy instead of
# re-downloading per-experiment. The catalog string is the version pin.
#
# NOTE: taxonomy.txt is NOT shared. It's a FILTERED projection of the
# metadata, with per-run filters (mgnify_taxonomy_filter, mgnify_max_genomes)
# already applied — different runs would produce different files. It lives
# in MGNIFY_OUT (per-run) alongside the per-run species_representatives.txt.
# Mixing per-run and shared outputs in the same checkpoint also confuses
# Snakemake's "is this checkpoint done?" check, which is why this matters.
MGNIFY_CACHE_DIR       = config.get("mgnify_cache_dir",
                                    "resources/genome_databases/mgnify")
MGNIFY_CATALOG_SLUG    = MGNIFY_CATALOG.replace("/", "_")
MGNIFY_CATALOG_ROOT    = (os.path.join(MGNIFY_CACHE_DIR, MGNIFY_CATALOG_SLUG)
                          if MGNIFY_CATALOG_SLUG else "")
MGNIFY_CACHE_GENOMES   = os.path.join(MGNIFY_CATALOG_ROOT, "genomes")
MGNIFY_CACHE_METADATA  = os.path.join(MGNIFY_CATALOG_ROOT, "genomes-all_metadata.tsv")

# ==============================================================================
# Helper functions
# ==============================================================================

def get_mgnify_genome_list():
    """Read species representative accessions from the checkpoint output."""
    reps_file = checkpoints.parse_mgnify_metadata.get().output.representatives
    with open(reps_file) as f:
        return sorted([line.strip() for line in f if line.strip()])


def mgnify_genome_path(accession):
    """Resolve a genome accession to its cached FASTA path."""
    return os.path.join(MGNIFY_CACHE_GENOMES, f"{accession}.fna")


def mgnify_genome_url(wildcards):
    """Build the FTP URL for a species representative genome.

    MGnify shards genomes on the FTP server by the first 11 chars of the
    accession (`MGYG` + 7 digits). Versioned accessions like
    `MGYG000000761.1` (15 chars) and unversioned ones like `MGYG000000761`
    (13 chars) live under the same 11-char bucket — the previous `acc[:-2]`
    slice happened to produce the right bucket for unversioned 13-char
    accessions but produced a non-existent path for any versioned one,
    causing silent 404 → HTML body saved as `.fna`.
    """
    acc = wildcards.accession
    prefix = acc[:11]
    return (
        f"{MGNIFY_FTP_BASE}/{MGNIFY_CATALOG}/"
        f"species_catalogue/{prefix}/{acc}/genome/{acc}.fna"
    )

# ==============================================================================
# Rules
# ==============================================================================

# Metadata lands in the shared catalog cache so every run against the same
# catalog reuses one download.
rule download_mgnify_metadata:
    output:
        metadata = MGNIFY_CACHE_METADATA
    params:
        url = f"{MGNIFY_FTP_BASE}/{MGNIFY_CATALOG}/genomes-all_metadata.tsv"
    log:
        os.path.join(MGNIFY_OUT, "logs/download_mgnify_metadata.log")
    shell:
        """
        mkdir -p $(dirname {output.metadata})
        mkdir -p $(dirname {log})
        curl -L --retry 5 --retry-delay 10 -C - \
            -o {output.metadata} '{params.url}' 2>&1 | tee {log}
        """


# Both outputs are per-run: taxonomy.txt is a filtered projection (depends on
# mgnify_taxonomy_filter + mgnify_max_genomes); representatives is the same
# filtered set as accessions only. Keeping both per-run keeps Snakemake's
# checkpoint "done?" check honest and avoids cross-run output collisions.
checkpoint parse_mgnify_metadata:
    input:
        metadata = MGNIFY_CACHE_METADATA
    output:
        taxonomy        = os.path.join(MGNIFY_OUT, "taxonomy.txt"),
        representatives = os.path.join(MGNIFY_OUT, "species_representatives.txt")
    params:
        taxonomy_filter = config.get("mgnify_taxonomy_filter", False),
        max_genomes     = config.get("mgnify_max_genomes", 0)
    log:
        os.path.join(MGNIFY_OUT, "logs/parse_mgnify_metadata.log")
    script:
        "scripts/parse_mgnify_metadata.py"


# Genomes land in the shared cache — two runs requesting the same accession
# Snakemake-dedupe to a single download.
rule download_mgnify_genome:
    input:
        representatives = os.path.join(MGNIFY_OUT, "species_representatives.txt")
    output:
        genome = os.path.join(MGNIFY_CACHE_GENOMES, "{accession}.fna")
    params:
        url = mgnify_genome_url
    log:
        os.path.join(MGNIFY_OUT, "logs/download_{accession}.log")
    # Run in the bakta container (same image the sibling genome-processing rules
    # use) so curl is a modern, reproducible version — the host curl is 7.68,
    # which lacks --retry-all-errors and made this rule host-version-dependent.
    container:
        config["containers"]["bakta"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.genome})
        mkdir -p $(dirname {log})
        # -f makes curl exit non-zero on 4xx/5xx (default behaviour writes the
        # HTML error body to {output.genome} and FragGeneScan downstream
        # segfaults on it). -S ensures errors are reported even with -s.
        # --retry-all-errors also retries connection-level failures (resets,
        # partial transfers) that plain --retry skips; --connect-timeout bounds
        # hangs against a flaky EBI FTP endpoint. (Requires the container's
        # modern curl; the host's 7.68 lacks --retry-all-errors.)
        curl -fSL --retry 8 --retry-delay 10 --retry-all-errors --connect-timeout 30 \
            -o {output.genome} '{params.url}' 2> >(tee -a {log} >&2)
        # Validate the downloaded file is actually a FASTA — defence-in-depth
        # in case the URL ever serves a 200 with junk content.
        if ! head -c 1 {output.genome} | grep -q '^>'; then
            echo "ERROR: {output.genome} does not start with '>' — not a FASTA" \
                | tee -a {log} >&2
            rm -f {output.genome}
            exit 1
        fi
        """


# Per-run sentinel — signals that every representative this run needs is
# present in the shared cache. Lives under MGNIFY_OUT (per-run) because the
# set of required genomes is a function of the per-run filtered list.
rule mgnify_download_complete:
    input:
        genomes = lambda wildcards: expand(
            os.path.join(MGNIFY_CACHE_GENOMES, "{acc}.fna"),
            acc=get_mgnify_genome_list()
        ),
        taxonomy = os.path.join(MGNIFY_OUT, "taxonomy.txt")
    output:
        touch(os.path.join(MGNIFY_OUT, ".mgnify_download_complete"))
