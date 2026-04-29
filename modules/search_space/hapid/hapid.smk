import os
import glob

include: "../_shared/genome_cache.smk"

# ==============================================================================
# Paths
# ==============================================================================
EXPERIMENT_DIR  = config["experiment_dir"]
RUN_DIR         = config["run_dir"]
HAPID_DIR       = os.path.join(EXPERIMENT_DIR, "input/MAG_files")
# Three tiers of cacheability — see _shared/genome_cache.smk:
#   HAPID_FGS_DIR    — per-genome (FragGeneScan FAA depends only on a single genome)
#   HAPID_HMMER_DIR  — per-genome (HMMER tblout depends only on FAA + HMM profiles)
#   HAPID_SET_DIR    — per-genome-set (deduped marker DB + speclib + protein2genome dict)
#   HAPID_OUT_ROOT   — per-run (DIA-NN profiling + greedy selection, MS-data-dependent)
HAPID_FGS_DIR   = per_genome_cache_root("hapid_fgs")
HAPID_HMMER_DIR = per_genome_cache_root("hapid_hmmer")
HAPID_SET_DIR   = per_genome_set_cache_root("hapid")
HAPID_OUT_ROOT  = os.path.join(RUN_DIR, "database_resources/hapid")

# MGnify shared cache (see modules/genome_download/mgnify/mgnify.smk).
_MGNIFY_CACHE_DIR    = config.get("mgnify_cache_dir",
                                  "resources/genome_databases/mgnify")
_MGNIFY_CATALOG_SLUG = config.get("mgnify_catalog", "").replace("/", "_")
_MGNIFY_CATALOG_ROOT = (os.path.join(_MGNIFY_CACHE_DIR, _MGNIFY_CATALOG_SLUG)
                        if _MGNIFY_CATALOG_SLUG else "")

def _mgnify_genome_path(genome):
    return os.path.join(_MGNIFY_CATALOG_ROOT, "genomes", f"{genome}.fna")

# ==============================================================================
# Helper functions
# ==============================================================================

# Re-import MGnify's representatives file through a checkpoint defined in this
# module — Snakemake 9 scopes the checkpoint proxy per-module, so we can't
# reach checkpoints.parse_mgnify_metadata directly from here. The input is the
# upstream file as a static path; checkpoint-aware DAG re-evaluation
# propagates correctly through it.
if config.get("genome_download_source") == "mgnify":
    checkpoint hapid_import_mgnify_genomes:
        input:
            os.path.join(RUN_DIR, "genome_download/mgnify/species_representatives.txt")
        output:
            os.path.join(HAPID_OUT_ROOT, "all_hapid_genomes.txt")
        log:
            os.path.join(RUN_DIR, "logs/search_space/hapid/import_mgnify_genomes.log")
        shell:
            "mkdir -p $(dirname {log}) && cp {input} {output} 2> {log}"

def get_all_hapid_genomes():
    """Genome IDs — from MGnify representatives when source=mgnify, else from
    user-provided FASTAs in HAPID_DIR. The mgnify branch must NOT fall through
    to the local scan: doing so caused atcc_25922 to be queued for MGnify
    download (404 → corrupt cache file) when the checkpoint hadn't fired yet."""
    if config.get("genome_download_source") == "mgnify":
        list_file = checkpoints.hapid_import_mgnify_genomes.get().output[0]
        with open(list_file) as f:
            return sorted([line.strip() for line in f if line.strip()])
    genomes = []
    for ext in ("fa", "fna", "fasta"):
        for f in glob.glob(os.path.join(HAPID_DIR, f"*.{ext}")):
            genomes.append(os.path.splitext(os.path.basename(f))[0])
    return sorted(set(genomes))


def hapid_fasta_path(wildcards):
    """Path to genome FASTA for a given wildcard. MGnify-sourced genomes come
    from the shared cache; user-provided genomes come from HAPID_DIR."""
    if config.get("genome_download_source") == "mgnify":
        return _mgnify_genome_path(wildcards.genome)
    for ext in ("fa", "fna", "fasta"):
        p = os.path.join(HAPID_DIR, f"{wildcards.genome}.{ext}")
        if os.path.exists(p):
            return p
    return os.path.join(HAPID_DIR, f"{wildcards.genome}.fa")


def marker_gene_db_path(wildcards=None):
    """Deduplicated marker gene FASTA."""
    return os.path.join(HAPID_SET_DIR, "marker_gene_db.fasta")


def marker_gene_clstr_path(wildcards=None):
    """CD-HIT cluster file."""
    return os.path.join(HAPID_SET_DIR, "marker_gene_db.fasta.clstr")


def protein2genome_dic_path(wildcards=None):
    """protein2genome JSON path."""
    return os.path.join(HAPID_SET_DIR, "protein2genome_dic.json")


# ==============================================================================
# Stage 1 — ORF prediction + HMMER marker gene identification
# ==============================================================================

rule press_hapid_hmm_profiles:
    input:
        config["hapid_hmm_profiles"]
    output:
        touch(config["hapid_hmm_profiles"] + ".pressed")
    log:
        os.path.join(RUN_DIR, "logs/search_space/hapid/press_hmm_profiles.log")
    container:
        config["containers"]["fraggenescan_hmmer"]
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
        os.path.join(HAPID_FGS_DIR, "{genome}.faa")
    params:
        prefix = lambda wildcards, output: output[0].replace(".faa", "")
    log:
        os.path.join(RUN_DIR, "logs/search_space/hapid/fgs/{genome}.log")
    threads: 1
    container:
        config["containers"]["fraggenescan_hmmer"]
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
        faa     = os.path.join(HAPID_FGS_DIR, "{genome}.faa"),
        pressed = config["hapid_hmm_profiles"] + ".pressed"
    output:
        os.path.join(HAPID_HMMER_DIR, "{genome}_hmmer.txt")
    params:
        hmm_profiles = config["hapid_hmm_profiles"]
    log:
        os.path.join(RUN_DIR, "logs/search_space/hapid/hmmer/{genome}.log")
    threads: 4
    container:
        config["containers"]["fraggenescan_hmmer"]
    shell:
        # hmmscan stdout is the verbose human-readable per-query report — we
        # already capture the parseable form via --tblout, so the stdout report
        # is dead weight (≈3 MB/genome × 4744 genomes ≈ 14 GB on a full MGnify
        # catalog). Discard it; keep stderr in {log} for real errors.
        """
        mkdir -p $(dirname {output})
        mkdir -p $(dirname {log})
        hmmscan \
            -E 1e-10 \
            --tblout {output} \
            --cpu {threads} \
            {params.hmm_profiles} \
            {input.faa} \
            > /dev/null 2> {log}
        """


rule build_hapid_marker_gene_fasta:
    input:
        hmmer_files = lambda wildcards: [
            os.path.join(HAPID_HMMER_DIR, f"{g}_hmmer.txt")
            for g in get_all_hapid_genomes()
        ],
        faa_files = lambda wildcards: [
            os.path.join(HAPID_FGS_DIR, f"{g}.faa")
            for g in get_all_hapid_genomes()
        ]
    output:
        os.path.join(HAPID_SET_DIR, "all_marker_genes.fasta")
    log:
        os.path.join(RUN_DIR, "logs/search_space/hapid/build_marker_gene_fasta.log")
    container:
        config["containers"]["fraggenescan_hmmer"]
    script:
        "scripts/build_marker_gene_fasta.py"


rule deduplicate_marker_genes_with_cdhit:
    input:
        fasta = os.path.join(HAPID_SET_DIR, "all_marker_genes.fasta")
    output:
        fasta = os.path.join(HAPID_SET_DIR, "marker_gene_db.fasta"),
        clstr = os.path.join(HAPID_SET_DIR, "marker_gene_db.fasta.clstr")
    log:
        os.path.join(RUN_DIR, "logs/search_space/hapid/cdhit.log")
    threads: 8
    container:
        config["containers"]["fraggenescan_hmmer"]
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
        fasta = os.path.join(HAPID_SET_DIR, "marker_gene_db.fasta")
    output:
        os.path.join(HAPID_SET_DIR, "protein2genome_dic.json")
    log:
        os.path.join(RUN_DIR, "logs/search_space/hapid/protein2genome_dict.log")
    container:
        config["containers"]["fraggenescan_hmmer"]
    script:
        "scripts/build_protein_genome_dict.py"


# ==============================================================================
# Stage 2 — DIA-NN profiling search
# ==============================================================================

rule create_hapid_profiling_spectral_library:
    input:
        fasta = marker_gene_db_path,
        cfg   = config["diann_spectral_library_base_config"]
    output:
        os.path.join(HAPID_SET_DIR, "marker_gene.predicted.speclib")
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
            --cut "K*,R*" \
            --missed-cleavages 1 \
            --min-pep-len 7 \
            --max-pep-len 30 \
            --species-ids \
            > {log} 2>&1
        """


rule perform_hapid_profiling_search:
    input:
        raw_dir = os.path.join(EXPERIMENT_DIR, "input/ms_files"),
        speclib = (
            [os.path.join(HAPID_SET_DIR, "marker_gene.predicted.speclib")]
            if config.get("hapid_search_mode", "standard") == "standard"
            else []
        ),
        fasta   = marker_gene_db_path,
        cfg     = (
            "config/hapid_infinidia.cfg"
            if config.get("hapid_search_mode", "standard") == "infinidia"
            else config["diann_library_search_base_config"]
        )
    output:
        os.path.join(HAPID_OUT_ROOT, "marker_gene_profiling_report.parquet")
    params:
        out_prefix = lambda wildcards, output: output[0].replace(".parquet", ""),
        lib_flag = (
            f"--lib {os.path.join(HAPID_SET_DIR, 'marker_gene.predicted.speclib')}"
            if config.get("hapid_search_mode", "standard") == "standard"
            else ""
        )
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
            {params.lib_flag} \
            --cut "K*,R*" \
            --missed-cleavages 1 \
            --min-pep-len 7 \
            --max-pep-len 30 \
            --threads {threads} \
            > {log} 2>&1
        """


# ==============================================================================
# Stage 3 — Greedy genome selection (checkpoint)
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
        config["containers"]["fraggenescan_hmmer"]
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
        config["containers"]["fraggenescan_hmmer"]
    shell:
        """
        mkdir -p $(dirname {log})
        python {workflow.basedir}/modules/search_space/_shared/scripts/coverAllSpectra_greedy.py \
            {input} {output} \
            > {log} 2>&1
        """


# Apply the hapid_percent_spectra cutoff and emit a one-genome-per-line list.
# Lives in this module so it has in-scope access to the greedy-selection
# checkpoint output; downstream MAGs rules consume the resulting file purely
# as a static input path.
rule hapid_filter_selected_genomes:
    input:
        os.path.join(HAPID_OUT_ROOT, "hapid_greedy_selection.tsv")
    output:
        os.path.join(HAPID_OUT_ROOT, "selected_genomes.txt")
    params:
        pct = config.get("hapid_percent_spectra", 80)
    run:
        import pandas as pd
        df = pd.read_csv(input[0], sep="\t")
        above = df[df["cumulative_pct"] >= params.pct]
        cutoff = (above.index[0] + 1) if not above.empty else len(df)
        with open(output[0], "w") as fh:
            for g in df["genome"].tolist()[:cutoff]:
                fh.write(f"{g}\n")
