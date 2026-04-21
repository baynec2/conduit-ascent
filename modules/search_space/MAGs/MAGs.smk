import os
import glob

# Experiment specific directories
EXPERIMENT_DIR = config["experiment_dir"]
RUN_DIR = config["run_dir"]
MAG_DIR = os.path.join(EXPERIMENT_DIR,"input/MAG_files")
# Resource specific directories.
BAKTA_DIR = config["bakta_db_dir"]
# Database specific output
BAKTA_OUT_ROOT = os.path.join(RUN_DIR,"database_resources/bakta")
DB_OUT_ROOT = os.path.join(RUN_DIR,"database_resources")

# Files that should be included in Bakta database
REQUIRED_BAKTA_FILES = (
    "bakta.db",
    "version.json",
    "expert-protein-sequences.dmnd",
    "sorf.dmnd",
    "psc.dmnd",
    "rfam-go.tsv",
    "oric.fna",
    "orit.fna",
)

# Get MAG names (basenames without extension) for wildcards
def get_mag_list():
    mags = []
    for ext in ("fa", "fna", "fasta"):
        for f in glob.glob(os.path.join(MAG_DIR, f"*.{ext}")):
            mags.append(os.path.splitext(os.path.basename(f))[0])
    return sorted(list(set(mags)))

# Get full path to MAG file given a MAG name (wildcard)
def mag_fasta_path(wildcards):
    for ext in ("fa", "fna", "fasta"):
        candidate = os.path.join(MAG_DIR, f"{wildcards.mag}.{ext}")
        if os.path.exists(candidate):
            return candidate
    return os.path.join(MAG_DIR, f"{wildcards.mag}.fa")

rule check_mag_fastas:
    input:
        MAG_DIR
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

bakta_db_final = os.path.join(f"{BAKTA_DIR}-{config['bakta_db_type']}")

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
        mag_fa = mag_fasta_path
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

rule create_uniprot_style_database:
    input:
        bakta_dirs = lambda wildcards: [
            os.path.join(BAKTA_OUT_ROOT, mag) for mag in get_mag_list()
        ]
    output:
        fasta = os.path.join(DB_OUT_ROOT, "mag_database.fasta"),
        go = os.path.join(DB_OUT_ROOT, "go_annotations.txt"),
        kegg = os.path.join(DB_OUT_ROOT, "kegg_annotations.txt")
    params:
        mag_metadata = os.path.join(MAG_DIR, "MAG_metadata.txt")
    log:
        os.path.join(RUN_DIR,"logs/search_space/MAGs/create_uniprot_sytle_database.log")
    container:
        config["containers"]["bakta"]
    script:
        "scripts/MAG_uniprot_headers.py"

rule get_mag_taxonomy:
    input:
        # File containing ncbi organism ids
        mag_metadata = os.path.join(EXPERIMENT_DIR,"input/MAG_files/MAG_metadata.txt")
    output:
        taxonomy = os.path.join(RUN_DIR,"database_resources/mag_taxonomy.txt")
    log: os.path.join(RUN_DIR,"logs/search_space/MAGs/get_mag_taxonomy.log")
    container: config["containers"]["conduitr"]
    script:
      "scripts/get_mag_taxonomy.R"

# This will allow us to integrate 
rule append_additional_organisms_or_proteomes:
    input:
        # Modifying the mag database to also have uniprot information. 
        mag_fasta = os.path.join(DB_OUT_ROOT, "mag_database.fasta"),
        mag_taxonomy = os.path.join(DB_OUT_ROOT,"mag_taxonomy.txt")
    output:
        # Directory containing uniprot annotations
        uniprot_fasta_dir = directory(os.path.join(DB_OUT_ROOT,"uniprot_database")),
        # modified files with 
        fasta = os.path.join(DB_OUT_ROOT,"database.fasta"),
        taxonomy = os.path.join(DB_OUT_ROOT,"taxonomy.txt")
    log: os.path.join(RUN_DIR,"logs/search_space/MAGs/append_additional_data.log")
    container: config["containers"]["conduitr"]
    script:
        "scripts/append_additional_organisms_or_proteomes.R"

rule get_mag_annotations:
    input:
        bakta_dirs = lambda wildcards: [
            os.path.join(BAKTA_OUT_ROOT, mag) for mag in get_mag_list()
        ],
    output:
        mag_annotations = os.path.join(RUN_DIR, "database_resources/bakta/mag_annotations.txt")
    log:
        os.path.join(RUN_DIR,"logs/search_space/MAGs/get_mag_annotations.log")
    container:
        config["containers"]["conduitr"]
    script:
        "scripts/get_mag_annotations.R"
