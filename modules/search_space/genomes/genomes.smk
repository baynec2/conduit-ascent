import os
import glob
import pandas as pd

include: "../_shared/genome_cache.smk"

# Experiment specific directories
EXPERIMENT_DIR = config["experiment_dir"]
RUN_DIR = config["run_dir"]
MAG_DIR = os.path.join(EXPERIMENT_DIR,"input/MAG_files")
# Resource specific directories.
BAKTA_DIR = config["bakta_db_dir"]
# Bakta annotations are shared across runs when source == "mgnify" (genome IDs are
# globally unique). For local MAGs they fall back to per-run paths to avoid
# cross-experiment name collisions. See _shared/genome_cache.smk.
BAKTA_OUT_ROOT = per_genome_cache_root("bakta")
DB_OUT_ROOT = os.path.join(RUN_DIR,"database_resources")

# MGnify shared cache paths (see modules/genome_download/mgnify/mgnify.smk).
# Duplicated here rather than imported because Snakemake modules don't share
# Python helpers across snakefiles; the logic is three lines.
_MGNIFY_CACHE_DIR    = config.get("mgnify_cache_dir",
                                  "resources/genome_databases/mgnify")
_MGNIFY_CATALOG_SLUG = config.get("mgnify_catalog", "").replace("/", "_")
_MGNIFY_CATALOG_ROOT = (os.path.join(_MGNIFY_CACHE_DIR, _MGNIFY_CATALOG_SLUG)
                        if _MGNIFY_CATALOG_SLUG else "")

def _mgnify_genome_path(mag):
    return os.path.join(_MGNIFY_CATALOG_ROOT, "genomes", f"{mag}.fna")

def _mgnify_taxonomy_path():
    # Per-run, NOT shared — the file is a filtered projection by per-run
    # mgnify_taxonomy_filter / mgnify_max_genomes, so two runs with different
    # filters would clobber each other in a shared location.
    return os.path.join(RUN_DIR, "genome_download/mgnify/taxonomy.txt")

# Files that should be included in Bakta database.
# The protein-sequence-cluster diamond DB differs by DB type: the `full` DB
# ships psc.dmnd (UniRef90); the `light` DB drops that tier and ships only
# pscc.dmnd (UniRef50 centroids). Declaring the wrong one as a required output
# makes Snakemake flag download_bakta_resources as failed (missing output)
# even though `bakta_db download` completed successfully.
_BAKTA_PSC_FILE = "psc.dmnd" if config["bakta_db_type"] == "full" else "pscc.dmnd"
REQUIRED_BAKTA_FILES = (
    "bakta.db",
    "version.json",
    "expert-protein-sequences.dmnd",
    "sorf.dmnd",
    _BAKTA_PSC_FILE,
    "rfam-go.tsv",
    "oric.fna",
    "orit.fna",
)

# Static dispatch to the upstream-emitted genome list. Each upstream selector
# is responsible for writing a one-genome-per-line file at a known path; from
# here we only see file paths, never peer modules' checkpoint proxies (which
# Snakemake 9 scopes per-module).
def _selected_genomes_source():
    method = config.get("search_space_method")
    if method == "genome_peptidotyping":
        return os.path.join(RUN_DIR, "database_resources/genome_peptidotyping/detected_genomes.txt")
    if method == "hapid":
        return os.path.join(RUN_DIR, "database_resources/hapid/selected_genomes.txt")
    if config.get("genome_download_source") == "mgnify":
        return os.path.join(RUN_DIR, "genome_download/mgnify/species_representatives.txt")
    return None

# Local checkpoint so get_mag_list() stays inside this module's checkpoints
# proxy (Snakemake 9 scopes it per-module). Output is temp() so the checkpoint
# re-runs every invocation — otherwise lambdas calling .get() resolve to "<TBD>"
# when its output persists from a prior run while upstream regenerates input.
if _selected_genomes_source() is not None:
    checkpoint canonicalize_selected_genomes:
        input:
            _selected_genomes_source()
        output:
            temp(os.path.join(DB_OUT_ROOT, "selected_genomes.txt"))
        log:
            os.path.join(RUN_DIR, "logs/search_space/MAGs/canonicalize_selected_genomes.log")
        shell:
            "mkdir -p $(dirname {log}) && cp {input} {output} 2> {log}"

# Get MAG names (basenames without extension) for wildcards
def get_mag_list():
    if _selected_genomes_source() is not None:
        list_file = checkpoints.canonicalize_selected_genomes.get().output[0]
        with open(list_file) as fh:
            return sorted([line.strip() for line in fh if line.strip()])
    # Plain MAGs method: use all user-provided FASTAs in MAG_DIR.
    mags = []
    for ext in ("fa", "fna", "fasta"):
        for f in glob.glob(os.path.join(MAG_DIR, f"*.{ext}")):
            mags.append(os.path.splitext(os.path.basename(f))[0])
    return sorted(set(mags))

# Get full path to MAG file given a MAG name (wildcard). MGnify-sourced
# genomes live in the shared cache; user-provided MAGs live in MAG_DIR.
def mag_fasta_path(wildcards):
    if config.get("genome_download_source") == "mgnify":
        return _mgnify_genome_path(wildcards.mag)
    for ext in ("fa", "fna", "fasta"):
        candidate = os.path.join(MAG_DIR, f"{wildcards.mag}.{ext}")
        if os.path.exists(candidate):
            return candidate
    return os.path.join(MAG_DIR, f"{wildcards.mag}.fa")

rule check_mag_fastas:
    input:
        MAG_DIR,
        *([os.path.join(RUN_DIR, "genome_download/mgnify/.mgnify_download_complete")]
          if config.get("genome_download_source") == "mgnify" else [])
    output:
        touch(os.path.join(MAG_DIR, ".fastas_checked"))
    log:
        os.path.join(MAG_DIR, "logs/check_mag_fastas.log")
    container: config["containers"]["bakta"]
    shell:
        r"""
        mkdir -p $(dirname {log})
        echo "Checking MAG FASTA files in {input}" > {log} 2>&1

        # Snakemake input files are already expanded as a space-separated list
        files=({input})

        if [ ${{#files[@]}} -eq 0 ]; then
            echo "ERROR: No MAG FASTA files found in {input}" | tee -a {log}
            exit 1
        fi

        echo "MAG FASTA files found:" >> {log}
        printf "%s\n" "${{files[@]}}" >> {log}

        touch {output}
        """

bakta_db_final = os.path.join(BAKTA_DIR, f"db-{config['bakta_db_type']}")

# Download the bakta resources if they do not exist at the user specified resource path.
rule download_bakta_resources:
    output:
        bakta_db_files =expand(os.path.join(bakta_db_final, "{file}"), file=REQUIRED_BAKTA_FILES),
        amrfinder_db = directory(os.path.join(bakta_db_final, "amrfinderplus-db"))
    log:
        os.path.join(BAKTA_DIR, "logs/download_bakta_resources.log")
    container:
        config["containers"]["bakta"]
    params:
        bakta_db_dir = config["bakta_db_dir"],
        bakta_db_type = config["bakta_db_type"]
    shell:
        """
        mkdir -p $(dirname {log})
        echo "Starting Bakta DB download..." > {log}
        
        # Download the database
        bakta_db download --output {params.bakta_db_dir} --type {params.bakta_db_type} >> {log} 2>&1

        echo "Bakta DB download finished!" >> {log}
        """
# Bakta db to search: resolve to absolute path so it works when Snakemake is run
# from any directory (e.g. experiments/CB019) and so the container sees the same path.


rule annotate_mags_with_bakta:
    input:
        mags_ok = os.path.join(MAG_DIR, ".fastas_checked"),
        mag_fa = mag_fasta_path,
        # Depend on the bakta DB so Snakemake schedules download_bakta_resources
        # when it's missing (and skips it when the DB is already staged at
        # bakta_db_dir). Without this the download rule is an orphan — nothing
        # requests its output, so the DB is never fetched and bakta fails.
        bakta_db = expand(os.path.join(bakta_db_final, "{file}"), file=REQUIRED_BAKTA_FILES)
    output:
        directory(os.path.join(BAKTA_OUT_ROOT, "{mag}"))
    params:
        bakta_db_final = bakta_db_final
    log:
        os.path.join(BAKTA_OUT_ROOT, "logs/{mag}_bakta.log")
    threads: workflow.cores
    container:
        config["containers"]["bakta"]
    shell:
        r"""
        mkdir -p $(dirname {log})
        echo "Annotating {input.mag_fa}" > {log}
        
        # Ensure writable temp directory for tRNAscan-SE and other tools
        export TMPDIR=$(mktemp -d -p /tmp)
        export TEMP=$TMPDIR
        export TMP=$TMPDIR
        
        # Suppress matplotlib cache warnings (optional)
        export MPLCONFIGDIR=$TMPDIR/matplotlib
        mkdir -p $MPLCONFIGDIR
        
        echo "Using temp directory: $TMPDIR" >> {log}

        bakta \
            --db {params.bakta_db_final} \
            --threads {threads} \
            --output {output} \
            --prefix {wildcards.mag} \
            {input.mag_fa} >> {log} 2>&1
        
        BAKTA_EXIT=$?
        
        # Cleanup
        rm -rf $TMPDIR
        
        if [ $BAKTA_EXIT -ne 0 ]; then
            echo "Bakta failed with exit code $BAKTA_EXIT" >> {log}
            exit $BAKTA_EXIT
        fi

        echo "Finished {wildcards.mag}" >> {log}
        """

def _mag_taxonomy_input():
    """Taxonomy source: shared cache when mgnify, experiment-local otherwise."""
    if config.get("genome_download_source") == "mgnify":
        return _mgnify_taxonomy_path()
    return os.path.join(MAG_DIR, "taxonomy.txt")

def _parse_mag_taxonomy_inputs(wildcards):
    # Snakemake input functions are always called with a wildcards arg, even
    # when the rule has no wildcards — accept and ignore it.
    inputs = {"taxonomy": _mag_taxonomy_input()}
    # When an upstream selector (hapid / genome_peptidotyping / mgnify reps) is
    # in play, filter taxonomy.txt to the selected subset so the final taxonomy
    # mirrors what's actually in database.fasta. Pure-MAGs runs have no
    # selector and the user-provided taxonomy already matches MAG_files/.
    if _selected_genomes_source() is not None:
        inputs["selected_genomes"] = os.path.join(DB_OUT_ROOT, "selected_genomes.txt")
    return inputs


rule parse_mag_taxonomy:
    input:
        unpack(_parse_mag_taxonomy_inputs)
    output:
        # Intermediate: consumed by create_uniprot_style_database (needs the
        # `genome` column for per-MAG species lookup) and by
        # append_additional_organisms_or_proteomes (drops `genome`, writes
        # taxonomy.txt). Snakemake auto-deletes after both finish.
        temp(os.path.join(DB_OUT_ROOT, "mag_taxonomy.txt"))
    log:
        os.path.join(RUN_DIR, "logs/search_space/MAGs/parse_mag_taxonomy.log")
    container:
        config["containers"]["bakta"]
    script:
        "scripts/parse_mag_taxonomy.py"

rule create_uniprot_style_database:
    input:
        bakta_dirs = lambda wildcards: [
            os.path.join(BAKTA_OUT_ROOT, mag) for mag in get_mag_list()
        ],
        taxonomy = os.path.join(DB_OUT_ROOT, "mag_taxonomy.txt")
    output:
        # Intermediate: consumed by append_additional_organisms_or_proteomes,
        # which writes the final database.fasta. Auto-deleted after.
        fasta = temp(os.path.join(DB_OUT_ROOT, "mag_database.fasta"))
    params:
        taxonomy = os.path.join(DB_OUT_ROOT, "mag_taxonomy.txt")
    log:
        os.path.join(RUN_DIR,"logs/search_space/MAGs/create_uniprot_style_database.log")
    container:
        config["containers"]["bakta"]
    script:
        "scripts/MAG_uniprot_headers.py"

rule append_additional_organisms_or_proteomes:
    input:
        # Modifying the mag database to also have uniprot information.
        mag_fasta = os.path.join(DB_OUT_ROOT, "mag_database.fasta"),
        mag_taxonomy = os.path.join(DB_OUT_ROOT,"mag_taxonomy.txt")
    output:
        fasta = os.path.join(DB_OUT_ROOT,"database.fasta"),
        taxonomy = os.path.join(DB_OUT_ROOT,"taxonomy.txt")
    log: os.path.join(RUN_DIR,"logs/search_space/MAGs/append_additional_data.log")
    # conduitR::get_proteome_ids_from_organism_ids() sizes its worker pool as
    # future::availableCores() - 1; guarantee >=2 cores so it never resolves to 0
    # workers (which would error) on a single-core allocation.
    threads: min(8, workflow.cores)
    container: config["containers"]["conduitr"]
    script:
        "scripts/append_additional_organisms_or_proteomes.R"

rule get_mag_annotations:
    input:
        bakta_dirs = lambda wildcards: [
            os.path.join(BAKTA_OUT_ROOT, mag) for mag in get_mag_list()
        ],
    output:
        # Per-run: filtered to this experiment's selected genomes (get_mag_list).
        # Don't write under BAKTA_OUT_ROOT — that path may be a shared cache.
        mag_annotations = os.path.join(DB_OUT_ROOT, "mag_annotations.txt")
    log:
        os.path.join(RUN_DIR,"logs/search_space/MAGs/get_mag_annotations.log")
    container:
        config["containers"]["conduitr"]
    script:
        "scripts/get_mag_annotations.R"
