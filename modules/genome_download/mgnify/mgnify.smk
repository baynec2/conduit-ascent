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

# ==============================================================================
# Helper functions
# ==============================================================================

def get_mgnify_genome_list():
    """Read species representative accessions from the checkpoint output."""
    reps_file = checkpoints.parse_mgnify_metadata.get().output.representatives
    with open(reps_file) as f:
        return sorted([line.strip() for line in f if line.strip()])


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

rule download_mgnify_metadata:
    output:
        metadata = os.path.join(MGNIFY_OUT, "genomes-all_metadata.tsv")
    params:
        url = f"{MGNIFY_FTP_BASE}/{MGNIFY_CATALOG}/genomes-all_metadata.tsv"
    log:
        os.path.join(MGNIFY_OUT, "logs/download_mgnify_metadata.log")
    shell:
        """
        mkdir -p $(dirname {output.metadata})
        curl -L --retry 5 --retry-delay 10 -C - \
            -o {output.metadata} '{params.url}' 2>&1 | tee {log}
        """


checkpoint parse_mgnify_metadata:
    input:
        metadata = os.path.join(MGNIFY_OUT, "genomes-all_metadata.tsv")
    output:
        taxonomy        = os.path.join(MAG_DIR, "taxonomy.txt"),
        representatives = os.path.join(MGNIFY_OUT, "species_representatives.txt")
    params:
        taxonomy_filter = config.get("mgnify_taxonomy_filter", False),
        max_genomes     = config.get("mgnify_max_genomes", 0)
    log:
        os.path.join(MGNIFY_OUT, "logs/parse_mgnify_metadata.log")
    script:
        "scripts/parse_mgnify_metadata.py"


rule download_mgnify_genome:
    input:
        representatives = os.path.join(MGNIFY_OUT, "species_representatives.txt")
    output:
        genome = os.path.join(MAG_DIR, "{accession}.fna")
    params:
        url = mgnify_genome_url
    log:
        os.path.join(MGNIFY_OUT, "logs/download_{accession}.log")
    shell:
        """
        mkdir -p $(dirname {output.genome})
        curl -L --retry 5 --retry-delay 10 \
            -o {output.genome} '{params.url}' 2>&1 | tee {log}
        """


rule mgnify_download_complete:
    input:
        genomes = lambda wildcards: expand(
            os.path.join(MAG_DIR, "{acc}.fna"),
            acc=get_mgnify_genome_list()
        ),
        taxonomy = os.path.join(MAG_DIR, "taxonomy.txt")
    output:
        touch(os.path.join(MAG_DIR, ".mgnify_download_complete"))
