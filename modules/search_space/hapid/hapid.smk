import os
import glob
import sys

include: "../_shared/genome_cache.smk"

# ==============================================================================
# Paths
# ==============================================================================
EXPERIMENT_DIR  = config["experiment_dir"]
RUN_DIR         = config["run_dir"]
# User-provided genome FASTAs (hapid selects a subset of these). Prefer
# input/genome_files/; fall back to the legacy input/MAG_files/ name (deprecated).
HAPID_DIR = os.path.join(EXPERIMENT_DIR, "input/genome_files")
if not os.path.isdir(HAPID_DIR) and os.path.isdir(os.path.join(EXPERIMENT_DIR, "input/MAG_files")):
    HAPID_DIR = os.path.join(EXPERIMENT_DIR, "input/MAG_files")

sys.path.insert(0, os.path.join(workflow.basedir, "modules", "_shared"))
from diann_staging import (
    list_raw_files,
    list_samples,
    raw_path_for_sample,
    stage3_symlink_commands,
)

RAW_FILEPATHS = list_raw_files(EXPERIMENT_DIR)
SAMPLES = list_samples(EXPERIMENT_DIR)
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
        # ancient(): .fastas_checked is a per-run sentinel (re-touched every run,
        # transitively via .mgnify_download_complete) but this FAA lands in the
        # SHARED per-genome cache. Without ancient() its fresh mtime re-runs
        # FragGeneScan for all genomes every run. It's purely an ordering guard
        # (ensures the genome FASTAs exist); the FAA content depends only on
        # genome_fa, so ignoring its timestamp is safe.
        fastas_ok = ancient(os.path.join(HAPID_DIR, ".fastas_checked"))
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
        # FragGeneScan reads training files relative to its install dir, so we
        # must `cd` into it first. Resolve all paths to absolute form *before*
        # the cd so they survive (works whether {output} is a per-run relative
        # path or an absolute mgnify-cache path).
        """
        INPUT_FA=$(realpath -m {input.genome_fa})
        OUTPUT=$(realpath -m {output})
        PREFIX=$(realpath -m {params.prefix})
        LOG=$(realpath -m {log})
        mkdir -p "$(dirname "$OUTPUT")"
        mkdir -p "$(dirname "$LOG")"
        cd "$(dirname "$(which FragGeneScan)")"
        FragGeneScan \
            -s "$INPUT_FA" \
            -o "$PREFIX" \
            -w 0 \
            -t complete \
            -p {threads} \
            > "$LOG" 2>&1
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
        # cd + basename for --out-lib so dots in the parent dir (catalog slug) don't get parsed as extensions.
        """
        INPUT_FA=$(realpath -m {input.fasta})
        CFG=$(realpath -m {input.cfg})
        LOG=$(realpath -m {log})
        OUT_DIR=$(dirname $(realpath -m {output}))
        OUT_BASE=$(basename {params.out_prefix})
        mkdir -p "$OUT_DIR"
        mkdir -p "$(dirname "$LOG")"
        cd "$OUT_DIR"
        {config[diann_cmd]} \
            --cfg "$CFG" \
            --fasta "$INPUT_FA" \
            --threads {threads} \
            --out-lib "$OUT_BASE" \
            --cut "K*,R*" \
            --missed-cleavages 1 \
            --min-pep-len 7 \
            --max-pep-len 30 \
            --species-ids \
            > "$LOG" 2>&1
        """


GENOME_HAPID_QUANTS = os.path.join(HAPID_OUT_ROOT, "profiling_quant_files")

# HAPiID-style marker-gene profiling search.
# Standard: 3-stage split. InfinDIA: monolithic. See modules/diann/diann.smk for rationale.
if config.get("hapid_search_mode", "standard") == "standard":

    rule genome_hapid_build_empirical_lib:
        input:
            raw_dir = os.path.join(EXPERIMENT_DIR, "input/ms_files"),
            speclib = os.path.join(HAPID_SET_DIR, "marker_gene.predicted.speclib"),
            fasta   = marker_gene_db_path,
            cfg     = os.path.join(RUN_DIR,"config/diann_library_search_base.cfg")
        output:
            empirical_lib = os.path.join(HAPID_OUT_ROOT, "marker_gene_empirical.parquet")
        params:
            out_lib = os.path.join(HAPID_OUT_ROOT, "marker_gene_empirical"),
            tmpdir = os.path.join(HAPID_OUT_ROOT, "build_empirical_quant_files"),
        log:
            os.path.join(RUN_DIR, "logs/search_space/hapid/build_empirical_lib.log")
        threads: workflow.cores
        container:
            config["containers"]["diann"]
        shell:
            """
            mkdir -p $(dirname {output.empirical_lib})
            mkdir -p $(dirname {log})
            rm -rf {params.tmpdir} && mkdir -p {params.tmpdir}
            {config[diann_cmd]} \
                --cfg {input.cfg} \
                --fasta {input.fasta} \
                --dir {input.raw_dir} \
                --temp {params.tmpdir} \
                --lib {input.speclib} \
                --gen-spec-lib \
                --rt-profiling \
                --out-lib {params.out_lib} \
                --out {params.tmpdir}/report \
                --cut "K*,R*" \
                --missed-cleavages 1 \
                --min-pep-len 7 \
                --max-pep-len 30 \
                --threads {threads} \
                > {log} 2>&1
            """

    rule genome_hapid_search_one_raw:
        input:
            empirical_lib = os.path.join(HAPID_OUT_ROOT, "marker_gene_empirical.parquet"),
            fasta = marker_gene_db_path,
            cfg = os.path.join(RUN_DIR,"config/diann_library_search_base.cfg"),
            raw = lambda w: raw_path_for_sample(EXPERIMENT_DIR, w.sample)
        output:
            quant = os.path.join(GENOME_HAPID_QUANTS, "{sample}.quant")
        params:
            tmpdir = lambda w: os.path.join(HAPID_OUT_ROOT, "profiling_quant_tmp", w.sample)
        log:
            os.path.join(RUN_DIR, "logs/search_space/hapid/search_one_raw.{sample}.log")
        threads: min(8, workflow.cores)
        container:
            config["containers"]["diann"]
        shell:
            """
            mkdir -p $(dirname {log})
            rm -rf {params.tmpdir} && mkdir -p {params.tmpdir} $(dirname {output.quant})
            {config[diann_cmd]} \
                --cfg {input.cfg} \
                --f {input.raw} \
                --lib {input.empirical_lib} \
                --fasta {input.fasta} \
                --temp {params.tmpdir} \
                --out {params.tmpdir}/per_run_report \
                --cut "K*,R*" \
                --missed-cleavages 1 \
                --min-pep-len 7 \
                --max-pep-len 30 \
                --threads {threads} \
                > {log} 2>&1
            mv {params.tmpdir}/*.quant {output.quant}
            rm -rf {params.tmpdir}
            """

    rule genome_hapid_combine:
        input:
            quants = expand(
                os.path.join(GENOME_HAPID_QUANTS, "{sample}.quant"),
                sample=SAMPLES
            ),
            empirical_lib = os.path.join(HAPID_OUT_ROOT, "marker_gene_empirical.parquet"),
            fasta = marker_gene_db_path,
            cfg = os.path.join(RUN_DIR,"config/diann_library_search_base.cfg"),
            raw_dir = os.path.join(EXPERIMENT_DIR, "input/ms_files")
        output:
            os.path.join(HAPID_OUT_ROOT, "marker_gene_profiling_report.parquet")
        params:
            tmpdir = os.path.join(HAPID_OUT_ROOT, "profiling_combine_tmp"),
            out_prefix = lambda wildcards, output: output[0].replace(".parquet", ""),
            symlink_cmds = stage3_symlink_commands(
                os.path.join(HAPID_OUT_ROOT, "profiling_combine_tmp"),
                RAW_FILEPATHS,
                GENOME_HAPID_QUANTS
            )
        log:
            os.path.join(RUN_DIR, "logs/search_space/hapid/combine.log")
        threads: workflow.cores
        container:
            config["containers"]["diann"]
        shell:
            """
            mkdir -p $(dirname {output[0]})
            mkdir -p $(dirname {log})
            rm -rf {params.tmpdir} && mkdir -p {params.tmpdir}
            {params.symlink_cmds}
            {config[diann_cmd]} \
                --cfg {input.cfg} \
                --dir {input.raw_dir} \
                --lib {input.empirical_lib} \
                --fasta {input.fasta} \
                --temp {params.tmpdir} \
                --use-quant \
                --out {params.out_prefix} \
                --cut "K*,R*" \
                --missed-cleavages 1 \
                --min-pep-len 7 \
                --max-pep-len 30 \
                --threads {threads} \
                > {log} 2>&1
            rm -rf {params.tmpdir}
            """

else:  # infinidia — monolithic

    rule genome_hapid_monolithic:
        input:
            raw_dir = os.path.join(EXPERIMENT_DIR, "input/ms_files"),
            fasta   = marker_gene_db_path,
            cfg     = os.path.join(RUN_DIR,"config/hapid_infinidia.cfg")
        output:
            os.path.join(HAPID_OUT_ROOT, "marker_gene_profiling_report.parquet")
        params:
            out_prefix = lambda wildcards, output: output[0].replace(".parquet", ""),
            tmpdir = os.path.join(HAPID_OUT_ROOT, "monolithic_quant_files"),
        log:
            os.path.join(RUN_DIR, "logs/search_space/hapid/monolithic.log")
        threads: workflow.cores
        container:
            config["containers"]["diann"]
        shell:
            """
            mkdir -p $(dirname {output[0]})
            mkdir -p $(dirname {log})
            rm -rf {params.tmpdir} && mkdir -p {params.tmpdir}
            {config[diann_cmd]} \
                --cfg {input.cfg} \
                --fasta {input.fasta} \
                --dir {input.raw_dir} \
                --temp {params.tmpdir} \
                --pre-search --pre-filter \
                --gen-spec-lib \
                --rt-profiling \
                --out {params.out_prefix} \
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
        os.path.join(HAPID_OUT_ROOT, "hapid_greedy_selection_raw.tsv")
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


# Enrich the greedy TSV with a `taxon_name` column. For MGnify-sourced
# genomes the per-run taxonomy.txt already maps each accession to its GTDB
# species; for user-supplied genomes no taxonomy source exists so taxon_name
# is left blank.
def _hapid_annotate_inputs(wildcards=None):
    inputs = {"raw": os.path.join(HAPID_OUT_ROOT, "hapid_greedy_selection_raw.tsv")}
    if config.get("genome_download_source") == "mgnify":
        inputs["taxonomy"] = os.path.join(RUN_DIR, "genome_download/mgnify/taxonomy.txt")
    return inputs


rule annotate_hapid_greedy_selection:
    input:
        unpack(_hapid_annotate_inputs)
    output:
        os.path.join(HAPID_OUT_ROOT, "hapid_greedy_selection.tsv")
    log:
        os.path.join(RUN_DIR, "logs/search_space/hapid/annotate_greedy_selection.log")
    run:
        import pandas as pd
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        raw = pd.read_csv(input.raw, sep="\t", dtype={"genome": str})
        taxonomy_path = getattr(input, "taxonomy", None)
        if taxonomy_path:
            # taxonomy.txt columns: genome, domain, kingdom, phylum, class, order, family, genus, species.
            # GTDB sometimes resolves only to genus (or higher) — fall back through
            # ranks so the user always sees the most-specific available name.
            rank_cols = ["species", "genus", "family", "order", "class", "phylum", "kingdom", "domain"]
            tax = pd.read_csv(taxonomy_path, sep="\t", dtype=str)
            tax["taxon_name"] = tax[rank_cols].bfill(axis=1).iloc[:, 0]
            merged = raw.merge(tax[["genome", "taxon_name"]], on="genome", how="left")
        else:
            merged = raw.copy()
            merged["taxon_name"] = ""
        merged = merged[["genome", "taxon_name", "nSpectraCovered", "cumulative_pct"]]
        merged.to_csv(output[0], sep="\t", index=False)
        with open(log[0], "w") as lh:
            n_missing = merged["taxon_name"].isna().sum() + (merged["taxon_name"] == "").sum()
            lh.write(
                f"Annotated {len(merged)} genomes with taxon_name "
                f"(source={'mgnify taxonomy.txt' if taxonomy_path else 'none'}, "
                f"{n_missing} blank).\n"
            )


# Apply the hapid_percent_spectra cutoff and emit a one-genome-per-line list.
# Lives in this module so it has in-scope access to the greedy-selection
# checkpoint output; downstream genomes rules consume the resulting file purely
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
