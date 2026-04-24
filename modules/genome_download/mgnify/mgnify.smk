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

# Globalized cache: genomes, metadata, and taxonomy live under a per-catalog
# directory keyed by the catalog string (e.g., "human-gut/v2.0.2" →
# "human-gut_v2.0.2") so multiple experiments / runs using the same catalog
# share one copy instead of re-downloading per-experiment. The catalog string
# is the version pin; no extra versioning is needed.
MGNIFY_CACHE_DIR       = config.get("mgnify_cache_dir",
                                    "resources/genome_databases/mgnify")
MGNIFY_CATALOG_SLUG    = MGNIFY_CATALOG.replace("/", "_")
MGNIFY_CATALOG_ROOT    = (os.path.join(MGNIFY_CACHE_DIR, MGNIFY_CATALOG_SLUG)
                          if MGNIFY_CATALOG_SLUG else "")
MGNIFY_CACHE_GENOMES   = os.path.join(MGNIFY_CATALOG_ROOT, "genomes")
MGNIFY_CACHE_METADATA  = os.path.join(MGNIFY_CATALOG_ROOT, "genomes-all_metadata.tsv")
MGNIFY_CACHE_TAXONOMY  = os.path.join(MGNIFY_CATALOG_ROOT, "taxonomy.txt")

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
    """Build the FTP URL for a species representative genome."""
    acc = wildcards.accession
    prefix = acc[:-2]
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


# Taxonomy is derived from the shared metadata and lives in the shared cache.
# species_representatives.txt stays per-run because mgnify_taxonomy_filter /
# mgnify_max_genomes may differ between runs.
checkpoint parse_mgnify_metadata:
    input:
        metadata = MGNIFY_CACHE_METADATA
    output:
        taxonomy        = MGNIFY_CACHE_TAXONOMY,
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
    shell:
        """
        mkdir -p $(dirname {output.genome})
        mkdir -p $(dirname {log})
        curl -L --retry 5 --retry-delay 10 \
            -o {output.genome} '{params.url}' 2>&1 | tee {log}
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
        taxonomy = MGNIFY_CACHE_TAXONOMY
    output:
        touch(os.path.join(MGNIFY_OUT, ".mgnify_download_complete"))
