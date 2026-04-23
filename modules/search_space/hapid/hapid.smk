import os
import glob
import pandas as pd

# ==============================================================================
# Paths and mode flags
# ==============================================================================
EXPERIMENT_DIR  = config["experiment_dir"]
RUN_DIR         = config["run_dir"]
HAPID_DIR       = os.path.join(EXPERIMENT_DIR, "input/MAG_files")
BAKTA_DIR       = config["bakta_db_dir"]
BAKTA_OUT_ROOT  = os.path.join(RUN_DIR, "database_resources/bakta")
DB_OUT_ROOT     = os.path.join(RUN_DIR, "database_resources")
HAPID_OUT_ROOT  = os.path.join(DB_OUT_ROOT, "hapid")
bakta_db_final  = os.path.join(BAKTA_DIR, f"db-{config['bakta_db_type']}")

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

# ==============================================================================
# Helper functions
# ==============================================================================

def get_all_hapid_genomes():
    """Genome IDs from user-provided genome FASTAs (genome mode only)."""
    if config.get("genome_download_source") == "mgnify":
        reps_file = checkpoints.parse_mgnify_metadata.get().output.representatives
        with open(reps_file) as f:
            return sorted([line.strip() for line in f if line.strip()])
    genomes = []
    for ext in ("fa", "fna", "fasta"):
        for f in glob.glob(os.path.join(HAPID_DIR, f"*.{ext}")):
            genomes.append(os.path.splitext(os.path.basename(f))[0])
    return sorted(set(genomes))


def hapid_fasta_path(wildcards):
    """Path to user-provided genome FASTA for a given wildcard."""
    for ext in ("fa", "fna", "fasta"):
        p = os.path.join(HAPID_DIR, f"{wildcards.genome}.{ext}")
        if os.path.exists(p):
            return p
    return os.path.join(HAPID_DIR, f"{wildcards.genome}.fa")


def get_selected_genomes(wildcards):
    """Return genome IDs covering hapid_percent_spectra% of profiling spectra."""
    chk = checkpoints.run_greedy_genome_selection.get().output[0]
    df  = pd.read_csv(chk, sep="\t")
    pct = config.get("hapid_percent_spectra", 80)
    above = df[df["cumulative_pct"] >= pct]
    cutoff = (above.index[0] + 1) if not above.empty else len(df)
    return df["genome"].tolist()[:cutoff]


def marker_gene_db_path(wildcards=None):
    """Deduplicated marker gene FASTA."""
    return os.path.join(HAPID_OUT_ROOT, "marker_gene_db.fasta")


def marker_gene_clstr_path(wildcards=None):
    """CD-HIT cluster file."""
    return os.path.join(HAPID_OUT_ROOT, "marker_gene_db.fasta.clstr")


def protein2genome_dic_path(wildcards=None):
    """protein2genome JSON path."""
    return os.path.join(HAPID_OUT_ROOT, "protein2genome_dic.json")


def hapid_database_fasta_path(wildcards=None):
    """Path of the hapid protein database before append step."""
    return os.path.join(DB_OUT_ROOT, "hapid_database.fasta")


# ==============================================================================
# Stage 1 — Genome mode: validation + FGS + scatter HMMER
# ==============================================================================

rule check_hapid_fastas:
    input:
        HAPID_DIR,
        *([os.path.join(HAPID_DIR, ".mgnify_download_complete")]
          if config.get("genome_download_source") == "mgnify" else [])
    output:
        touch(os.path.join(HAPID_DIR, ".fastas_checked"))
    log:
        os.path.join(HAPID_DIR, "logs/check_hapid_fastas.log")
    container:
        config["containers"]["bakta"]
    shell:
        r"""
        mkdir -p $(dirname {log})
        echo "Checking HAPiID genome FASTA files in {input}" > {log} 2>&1

        shopt -s nullglob
        files=({input}/*.fa {input}/*.fna {input}/*.fasta)

        if [ ${{#files[@]}} -eq 0 ]; then
            echo "ERROR: No genome FASTA files found in {input}" | tee -a {log}
            exit 1
        fi

        echo "Genome FASTA files found:" >> {log}
        printf "%s\n" "${{files[@]}}" >> {log}

        touch {output}
        """


rule press_hapid_hmm_profiles:
    input:
        config["hapid_hmm_profiles"]
    output:
        touch(config["hapid_hmm_profiles"] + ".pressed")
    log:
        os.path.join(RUN_DIR, "logs/search_space/hapid/press_hmm_profiles.log")
    container:
        config["hapid_fgs_hmmer_container"]
    shell:
        """
        mkdir -p $(dirname {log})
        hmmpress -f {input} > {log} 2>&1
        touch {output}
        """


rule predict_orfs_with_fraggenescan:
    input:
        genome_fa = hapid_fasta_path,
        fastas_ok = os.path.join(HAPID_DIR, ".fastas_checked")
    output:
        os.path.join(HAPID_OUT_ROOT, "fgs/{genome}.faa")
    params:
        prefix = lambda wildcards, output: output[0].replace(".faa", "")
    log:
        os.path.join(RUN_DIR, "logs/search_space/hapid/fgs/{genome}.log")
    threads: 1
    container:
        config["hapid_fgs_hmmer_container"]
    shell:
        """
        WORKDIR=$(pwd)
        mkdir -p $WORKDIR/$(dirname {output})
        mkdir -p $WORKDIR/$(dirname {log})
        cd $(dirname $(which FragGeneScan))
        FragGeneScan \
            -s $WORKDIR/{input.genome_fa} \
            -o $WORKDIR/{params.prefix} \
            -w 0 \
            -t complete \
            -p {threads} \
            > $WORKDIR/{log} 2>&1
        """


rule identify_marker_genes_with_hmmer:
    input:
        faa     = os.path.join(HAPID_OUT_ROOT, "fgs/{genome}.faa"),
        pressed = config["hapid_hmm_profiles"] + ".pressed"
    output:
        os.path.join(HAPID_OUT_ROOT, "hmmer/{genome}_hmmer.txt")
    params:
        hmm_profiles = config["hapid_hmm_profiles"]
    log:
        os.path.join(RUN_DIR, "logs/search_space/hapid/hmmer/{genome}.log")
    threads: 4
    container:
        config["hapid_fgs_hmmer_container"]
    shell:
        """
        mkdir -p $(dirname {output})
        mkdir -p $(dirname {log})
        hmmscan \
            -E 1e-10 \
            --tblout {output} \
            --cpu {threads} \
            {params.hmm_profiles} \
            {input.faa} \
            > {log} 2>&1
        """


rule build_hapid_marker_gene_fasta:
    input:
        hmmer_files = lambda wildcards: [
            os.path.join(HAPID_OUT_ROOT, f"hmmer/{g}_hmmer.txt")
            for g in get_all_hapid_genomes()
        ],
        faa_files = lambda wildcards: [
            os.path.join(HAPID_OUT_ROOT, f"fgs/{g}.faa")
            for g in get_all_hapid_genomes()
        ]
    output:
        os.path.join(HAPID_OUT_ROOT, "all_marker_genes.fasta")
    log:
        os.path.join(RUN_DIR, "logs/search_space/hapid/build_marker_gene_fasta.log")
    container:
        config["hapid_fgs_hmmer_container"]
    script:
        "scripts/build_marker_gene_fasta.py"


rule deduplicate_marker_genes_with_cdhit:
    input:
        fasta = os.path.join(HAPID_OUT_ROOT, "all_marker_genes.fasta")
    output:
        fasta = os.path.join(HAPID_OUT_ROOT, "marker_gene_db.fasta"),
        clstr = os.path.join(HAPID_OUT_ROOT, "marker_gene_db.fasta.clstr")
    log:
        os.path.join(RUN_DIR, "logs/search_space/hapid/cdhit.log")
    threads: 8
    container:
        config["hapid_fgs_hmmer_container"]
    shell:
        """
        mkdir -p $(dirname {log})
        cd-hit \
            -i {input.fasta} \
            -o {output.fasta} \
            -c 1.0 \
            -n 5 \
            -T {threads} \
            -d 0 \
            -M 30000 \
            > {log} 2>&1
        """


rule build_protein_genome_dict:
    input:
        fasta = os.path.join(HAPID_OUT_ROOT, "marker_gene_db.fasta")
    output:
        os.path.join(HAPID_OUT_ROOT, "protein2genome_dic.json")
    log:
        os.path.join(RUN_DIR, "logs/search_space/hapid/protein2genome_dict.log")
    container:
        config["hapid_fgs_hmmer_container"]
    script:
        "scripts/build_protein_genome_dict.py"


# ==============================================================================
# Stage 2 — DIA-NN profiling search (shared by all modes)
# ==============================================================================

rule create_hapid_profiling_spectral_library:
    input:
        fasta = marker_gene_db_path,
        cfg   = "config/hapid_profiling_diann_spectral_library.cfg"
    output:
        os.path.join(HAPID_OUT_ROOT, "marker_gene.predicted.speclib")
    params:
        out_prefix = lambda wildcards, output: output[0].replace(".predicted.speclib", "")
    log:
        os.path.join(RUN_DIR, "logs/search_space/hapid/profiling_speclib.log")
    threads: workflow.cores
    container:
        config["containers"]["diann"]
    shell:
        """
        mkdir -p $(dirname {output})
        mkdir -p $(dirname {log})
        diann \
            --cfg {input.cfg} \
            --fasta {input.fasta} \
            --threads {threads} \
            --out-lib {params.out_prefix} \
            > {log} 2>&1
        """


rule perform_hapid_profiling_search:
    input:
        raw_dir = os.path.join(EXPERIMENT_DIR, "input/raw_files"),
        speclib = os.path.join(HAPID_OUT_ROOT, "marker_gene.predicted.speclib"),
        fasta   = marker_gene_db_path,
        cfg     = "config/hapid_profiling_diann.cfg"
    output:
        os.path.join(HAPID_OUT_ROOT, "marker_gene_profiling_report.parquet")
    params:
        out_prefix = lambda wildcards, output: output[0].replace(".parquet", "")
    log:
        os.path.join(RUN_DIR, "logs/search_space/hapid/profiling_search.log")
    threads: workflow.cores
    container:
        config["containers"]["diann"]
    shell:
        """
        mkdir -p $(dirname {output})
        mkdir -p $(dirname {log})
        diann \
            --cfg {input.cfg} \
            --fasta {input.fasta} \
            --out {params.out_prefix} \
            --dir {input.raw_dir} \
            --lib {input.speclib} \
            --threads {threads} \
            > {log} 2>&1
        """


# ==============================================================================
# Stage 3 — Greedy genome selection (checkpoint) — shared by all modes
# ==============================================================================

rule build_genome_spectrum_mapping:
    input:
        parquet        = os.path.join(HAPID_OUT_ROOT, "marker_gene_profiling_report.parquet"),
        protein2genome = protein2genome_dic_path,
        clstr          = marker_gene_clstr_path
    output:
        os.path.join(HAPID_OUT_ROOT, "genome2spectrum_dic.json")
    log:
        os.path.join(RUN_DIR, "logs/search_space/hapid/build_genome_spectrum_mapping.log")
    container:
        config["hapid_fgs_hmmer_container"]
    script:
        "scripts/build_genome_spectrum_mapping.py"


checkpoint run_greedy_genome_selection:
    input:
        os.path.join(HAPID_OUT_ROOT, "genome2spectrum_dic.json")
    output:
        os.path.join(HAPID_OUT_ROOT, "hapid_greedy_selection.tsv")
    log:
        os.path.join(RUN_DIR, "logs/search_space/hapid/greedy_selection.log")
    container:
        config["hapid_fgs_hmmer_container"]
    shell:
        """
        mkdir -p $(dirname {log})
        python {workflow.basedir}/modules/search_space/hapid/scripts/coverAllSpectra_greedy.py \
            {input} {output} \
            > {log} 2>&1
        """


# ==============================================================================
# Stage 3b — Filter taxonomy to selected genomes
# ==============================================================================

rule parse_hapid_selected_taxonomy:
    input:
        taxonomy  = os.path.join(HAPID_DIR, "taxonomy.txt"),
        selection = os.path.join(HAPID_OUT_ROOT, "hapid_greedy_selection.tsv")
    output:
        os.path.join(DB_OUT_ROOT, "hapid_selected_taxonomy.txt")
    log:
        os.path.join(RUN_DIR, "logs/search_space/hapid/parse_hapid_selected_taxonomy.log")
    params:
        hapid_percent_spectra = config.get("hapid_percent_spectra", 80)
    container:
        config["hapid_fgs_hmmer_container"]
    script:
        "scripts/parse_hapid_selected_taxonomy.py"


# ==============================================================================
# Stage 4 — Bakta annotation on selected genomes
# ==============================================================================

rule download_bakta_resources:
    output:
        bakta_db_files = expand(
            os.path.join(bakta_db_final, "{file}"), file=REQUIRED_BAKTA_FILES
        ),
        amrfinder_db = directory(os.path.join(bakta_db_final, "amrfinderplus-db"))
    log:
        os.path.join(BAKTA_DIR, "logs/download_bakta_resources.log")
    container:
        config["containers"]["bakta"]
    params:
        bakta_db_dir  = config["bakta_db_dir"],
        bakta_db_type = config["bakta_db_type"]
    shell:
        """
        mkdir -p $(dirname {log})
        echo "Starting Bakta DB download..." > {log}
        bakta_db download --output {params.bakta_db_dir} --type {params.bakta_db_type} >> {log} 2>&1
        echo "Bakta DB download finished!" >> {log}
        """


rule annotate_selected_hapid_genomes_with_bakta:
    input:
        fastas_ok = os.path.join(HAPID_DIR, ".fastas_checked"),
        genome_fa = hapid_fasta_path,
        bakta_db  = expand(
            os.path.join(bakta_db_final, "{file}"), file=REQUIRED_BAKTA_FILES
        )
    output:
        directory(os.path.join(BAKTA_OUT_ROOT, "{genome}"))
    params:
        bakta_db_final = bakta_db_final
    log:
        os.path.join(BAKTA_OUT_ROOT, "logs/{genome}_bakta.log")
    threads: workflow.cores
    container:
        config["containers"]["bakta"]
    shell:
        r"""
        mkdir -p $(dirname {log})
        echo "Annotating {input.genome_fa}" > {log}

        export TMPDIR=$(mktemp -d -p /tmp)
        export TEMP=$TMPDIR
        export TMP=$TMPDIR
        export MPLCONFIGDIR=$TMPDIR/matplotlib
        mkdir -p $MPLCONFIGDIR

        bakta \
            --db {params.bakta_db_final} \
            --threads {threads} \
            --output {output} \
            --prefix {wildcards.genome} \
            {input.genome_fa} >> {log} 2>&1

        BAKTA_EXIT=$?
        rm -rf $TMPDIR

        if [ $BAKTA_EXIT -ne 0 ]; then
            echo "Bakta failed with exit code $BAKTA_EXIT" >> {log}
            exit $BAKTA_EXIT
        fi

        echo "Finished {wildcards.genome}" >> {log}
        """


rule create_hapid_uniprot_style_database:
    input:
        bakta_dirs     = lambda wildcards: [
            os.path.join(BAKTA_OUT_ROOT, g) for g in get_selected_genomes(wildcards)
        ],
        hapid_taxonomy = os.path.join(DB_OUT_ROOT, "hapid_selected_taxonomy.txt")
    output:
        fasta = os.path.join(DB_OUT_ROOT, "hapid_database.fasta"),
        go    = os.path.join(DB_OUT_ROOT, "go_annotations.txt"),
        kegg  = os.path.join(DB_OUT_ROOT, "kegg_annotations.txt")
    params:
        hapid_taxonomy = os.path.join(DB_OUT_ROOT, "hapid_selected_taxonomy.txt")
    log:
        os.path.join(RUN_DIR, "logs/search_space/hapid/create_uniprot_style_database.log")
    container:
        config["containers"]["bakta"]
    script:
        "scripts/hapid_uniprot_headers.py"


# ==============================================================================
# Stage 5 — Taxonomy, annotations, final database (shared, mode-aware inputs)
# ==============================================================================

rule parse_hapid_taxonomy:
    input:
        taxonomy = os.path.join(HAPID_DIR, "taxonomy.txt")
    output:
        os.path.join(DB_OUT_ROOT, "hapid_taxonomy.txt")
    log:
        os.path.join(RUN_DIR, "logs/search_space/hapid/parse_hapid_taxonomy.log")
    container:
        config["hapid_fgs_hmmer_container"]
    script:
        "scripts/parse_hapid_taxonomy.py"


rule plot_input_taxonomy:
    input:
        os.path.join(DB_OUT_ROOT, "hapid_taxonomy.txt")
    output:
        os.path.join(HAPID_OUT_ROOT, "input_taxonomy_plot.pdf")
    log:
        os.path.join(RUN_DIR, "logs/search_space/hapid/plot_input_taxonomy.log")
    container:
        config["containers"]["conduitr"]
    script:
        "../database_processing/scripts/plot_taxonomic_tree.R"


rule get_hapid_annotations:
    input:
        bakta_dirs = lambda wildcards: [
            os.path.join(BAKTA_OUT_ROOT, g) for g in get_selected_genomes(wildcards)
        ]
    output:
        mag_annotations = os.path.join(
            RUN_DIR, "database_resources/bakta/mag_annotations.txt"
        )
    log:
        os.path.join(RUN_DIR, "logs/search_space/hapid/get_hapid_annotations.log")
    container:
        config["hapid_fgs_hmmer_container"]
    script:
        "scripts/get_hapid_annotations.py"


rule append_hapid_additional_organisms_or_proteomes:
    input:
        mag_fasta    = hapid_database_fasta_path,
        mag_taxonomy = os.path.join(DB_OUT_ROOT, "hapid_selected_taxonomy.txt")
    output:
        uniprot_fasta_dir = directory(os.path.join(DB_OUT_ROOT, "uniprot_database")),
        fasta             = os.path.join(DB_OUT_ROOT, "database.fasta"),
        taxonomy          = os.path.join(DB_OUT_ROOT, "taxonomy.txt")
    log:
        os.path.join(RUN_DIR, "logs/search_space/hapid/append_additional_data.log")
    container:
        config["containers"]["conduitr"]
    script:
        "../MAGs/scripts/append_additional_organisms_or_proteomes.R"
