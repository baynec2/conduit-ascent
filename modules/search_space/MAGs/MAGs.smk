import os
import glob

MAG_DIR        = "experiments/example/input/MAG_files"
BAKTA_DB       = "resources/bakta/db"
BAKTA_OUT_ROOT = "experiments/example/output/bakta"
DB_OUT_ROOT    = "experiments/example/output/database_resources"

def get_mag_list():
    mags = []
    for ext in ("fa", "fna", "fasta"):
        for f in glob.glob(os.path.join(MAG_DIR, f"*.{ext}")):
            mags.append(os.path.splitext(os.path.basename(f))[0])
    return sorted(list(set(mags)))

def mag_fasta_path(wildcards):
    for ext in ("fa", "fna", "fasta"):
        candidate = os.path.join(MAG_DIR, f"{wildcards.mag}.{ext}")
        if os.path.exists(candidate):
            return candidate
    return os.path.join(MAG_DIR, f"{wildcards.mag}.fa")

def get_bakta_required_files():
    required_full = [
        "taxonomy",
        "database.json",
        "sequences.fna",
        "proteins.faa"
    ]
    required_lite = [
        "taxonomy",
        "database.json"
    ]
    return required_full, required_lite


rule check_mag_fastas:
    input:
        MAG_DIR
    output:
        touch(os.path.join(MAG_DIR, ".fastas_checked"))
    log:
        os.path.join(MAG_DIR, "logs/check_mag_fastas.log")
    shell:
        r"""
        mkdir -p $(dirname {log})
        echo "Checking MAG FASTA files in {input}" > {log} 2>&1

        shopt -s nullglob

        files=({input}/*.fa {input}/*.fna {input}/*.fasta)

        # Snakemake-safe array length check
        if [ ${{#files[@]}} -eq 0 ]; then
            echo "ERROR: No MAG FASTA files found in {input}" | tee -a {log}
            exit 1
        fi

        echo "MAG FASTA files found:" >> {log}
        printf "%s\n" "${{files[@]}}" >> {log}

        touch {output}
        """


rule check_bakta_resources:
    input:
        database_dir = BAKTA_DB
    output:
        touch(os.path.join(BAKTA_DB, ".db_checked"))
    log:
        os.path.join(BAKTA_DB, "logs/check_bakta_resources.log")
    container:
        "docker://baynec2/bakta:alpha"
    shell:
        r"""
        mkdir -p $(dirname {log})
        echo "Checking Bakta v6 DB in {input.database_dir}" > {log}

        REQUIRED_FILES=(
            "bakta.db"
            "version.json"
            "expert-protein-sequences.dmnd"
            "sorf.dmnd"
            "psc.dmnd"
            "rfam-go.tsv"
            "oric.fna"
            "orit.fna"
        )

        REQUIRED_DIRS=(
            "amrfinderplus-db"
        )

        for f in "${{REQUIRED_FILES[@]}}"; do
            if [ ! -e "{input.database_dir}/$f" ]; then
                echo "Missing required Bakta v6 file: $f" >> {log}
                exit 1
            fi
        done

        for d in "${{REQUIRED_DIRS[@]}}"; do
            if [ ! -d "{input.database_dir}/$d" ]; then
                echo "Missing required Bakta v6 directory: $d" >> {log}
                exit 1
            fi
        done

        echo "Bakta v6 database OK." >> {log}
        touch {output}
        """


rule annotate_mags_with_bakta:
    input:
        db_ok   = os.path.join(BAKTA_DB, ".db_checked"),
        mags_ok = os.path.join(MAG_DIR, ".fastas_checked"),
        mag_fa  = mag_fasta_path,
    output:
        directory(os.path.join(BAKTA_OUT_ROOT, "{mag}"))
    params:
        bakta_db = BAKTA_DB
    log:
        os.path.join(BAKTA_OUT_ROOT, "logs/{mag}_bakta.log")
    threads: 4
    container:
        "docker://baynec2/bakta:alpha"
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
            --db {params.bakta_db} \
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
        fasta = os.path.join(DB_OUT_ROOT, "MAGS_uniprot.fasta"),
        go    = os.path.join(DB_OUT_ROOT, "detected_protein_resources/MAGS_go_annotations.txt"),
        kegg  = os.path.join(DB_OUT_ROOT, "detected_protein_resources/MAGS_kegg_annotations.txt")
    log:
        os.path.join(DB_OUT_ROOT, "logs/create_uniprot_headers.log")
    script:
        "modules/search_space/MAGs/scripts/MAG_uniprot_headers.py"
