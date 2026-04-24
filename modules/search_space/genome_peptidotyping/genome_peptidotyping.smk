import glob
import os

EXPERIMENT_DIR = config["experiment_dir"]
RUN_DIR = config["run_dir"]
MAG_DIR = os.path.join(EXPERIMENT_DIR, "input/MAG_files")
RAW_FILEPATHS = glob.glob(os.path.join(EXPERIMENT_DIR, "input/ms_files/*.raw"))

# Output directories for genome peptidotyping intermediate files
GP_RESOURCE_DIR = os.path.join(RUN_DIR, "database_resources/genome_peptidotyping")
GP_PRODIGAL_DIR = os.path.join(GP_RESOURCE_DIR, "prodigal")

# MGnify shared cache (see modules/genome_download/mgnify/mgnify.smk).
_MGNIFY_CACHE_DIR    = config.get("mgnify_cache_dir",
                                  "resources/genome_databases/mgnify")
_MGNIFY_CATALOG_SLUG = config.get("mgnify_catalog", "").replace("/", "_")
_MGNIFY_CATALOG_ROOT = (os.path.join(_MGNIFY_CACHE_DIR, _MGNIFY_CATALOG_SLUG)
                        if _MGNIFY_CATALOG_SLUG else "")

def _mgnify_genome_path(genome):
    return os.path.join(_MGNIFY_CATALOG_ROOT, "genomes", f"{genome}.fna")

def _mgnify_taxonomy_path():
    return os.path.join(_MGNIFY_CATALOG_ROOT, "taxonomy.txt")

def _taxonomy_input():
    """Taxonomy source: shared cache when mgnify, experiment-local otherwise."""
    if config.get("genome_download_source") == "mgnify":
        return _mgnify_taxonomy_path()
    return os.path.join(MAG_DIR, "taxonomy.txt")

# Get all genome names from MAG_files directory
def get_all_genome_names():
    if config.get("genome_download_source") == "mgnify":
        reps_file = checkpoints.parse_mgnify_metadata.get().output.representatives
        if os.path.exists(reps_file):
            with open(reps_file) as f:
                return sorted([line.strip() for line in f if line.strip()])
    genomes = []
    for ext in ("fa", "fna", "fasta"):
        for f in glob.glob(os.path.join(MAG_DIR, f"*.{ext}")):
            genomes.append(os.path.splitext(os.path.basename(f))[0])
    return sorted(list(set(genomes)))

# Get full path to genome FASTA given a genome name. MGnify-sourced genomes
# live in the shared cache; user-provided genomes live in MAG_DIR.
def genome_fasta_path(wildcards):
    if config.get("genome_download_source") == "mgnify":
        return _mgnify_genome_path(wildcards.genome)
    for ext in ("fa", "fna", "fasta"):
        candidate = os.path.join(MAG_DIR, f"{wildcards.genome}.{ext}")
        if os.path.exists(candidate):
            return candidate
    return os.path.join(MAG_DIR, f"{wildcards.genome}.fa")

# Rank config — same as peptidotyping
GENOME_PEPTIDOTYPING_RANK_CONFIG = {
    "species_strain": "species,strain",
    "genus": "genus",
    "family": "family",
}

################################################################################
# Phase A: Lightweight Gene Prediction with Prodigal
################################################################################

rule predict_orfs_with_prodigal:
    input:
        genome_fa = genome_fasta_path
    output:
        faa = os.path.join(GP_PRODIGAL_DIR, "{genome}.faa")
    log:
        os.path.join(GP_RESOURCE_DIR, "logs/prodigal/{genome}.log")
    container:
        config["containers"]["bakta"]
    shell:
        r"""
        mkdir -p $(dirname {output.faa})
        mkdir -p $(dirname {log})
        pyrodigal -i {input.genome_fa} -a {output.faa} -p meta >> {log} 2>&1
        """

################################################################################
# Phase B: Tryptic Digest + LCA Computation
################################################################################

rule tryptic_digest_genomes:
    input:
        faa_files = lambda wildcards: [
            os.path.join(GP_PRODIGAL_DIR, f"{genome}.faa")
            for genome in get_all_genome_names()
        ],
        taxonomy = _taxonomy_input()
    output:
        peptide_mapping = os.path.join(GP_RESOURCE_DIR, "peptide_genome_mapping.tsv.gz")
    params:
        prodigal_dir = GP_PRODIGAL_DIR
    log:
        os.path.join(GP_RESOURCE_DIR, "logs/tryptic_digest.log")
    container:
        config["containers"]["bakta"]
    script:
        "scripts/tryptic_digest.py"


rule compute_peptide_lca_and_build_dbs:
    input:
        peptide_mapping = os.path.join(GP_RESOURCE_DIR, "peptide_genome_mapping.tsv.gz"),
        taxonomy = _taxonomy_input()
    output:
        family_tsv      = os.path.join(GP_RESOURCE_DIR, "family_lca_filtered_peptides.tsv"),
        family_fasta    = os.path.join(GP_RESOURCE_DIR, "family_peptidotyping_db.fasta"),
        genus_tsv       = os.path.join(GP_RESOURCE_DIR, "genus_lca_filtered_peptides.tsv"),
        genus_fasta     = os.path.join(GP_RESOURCE_DIR, "genus_peptidotyping_db.fasta"),
        species_tsv     = os.path.join(GP_RESOURCE_DIR, "species_strain_lca_filtered_peptides.tsv"),
        species_fasta   = os.path.join(GP_RESOURCE_DIR, "species_strain_peptidotyping_db.fasta"),
        taxid_family_map = os.path.join(GP_RESOURCE_DIR, "taxid_to_family_genus.tsv")
    log:
        os.path.join(GP_RESOURCE_DIR, "logs/compute_peptide_lca.log")
    container:
        config["containers"]["bakta"]
    script:
        "scripts/compute_peptide_lca.py"


################################################################################
# Phase C: Build Effective Detection Rank Database
################################################################################
# Adapted from peptidotyping's build_effective_detection_rank_db.
# Simplified: no taxonkit needed — taxid_to_family_genus.tsv was built in Phase B.

rule build_genome_peptidotyping_effective_detection_rank_db:
    input:
        family_tsv       = os.path.join(GP_RESOURCE_DIR, "family_lca_filtered_peptides.tsv"),
        genus_tsv        = os.path.join(GP_RESOURCE_DIR, "genus_lca_filtered_peptides.tsv"),
        species_tsv      = os.path.join(GP_RESOURCE_DIR, "species_strain_lca_filtered_peptides.tsv"),
        taxid_family_map = os.path.join(GP_RESOURCE_DIR, "taxid_to_family_genus.tsv")
    output:
        first_pass_fasta = os.path.join(GP_RESOURCE_DIR, "effective_first_pass_database.fasta"),
        rank_mapping     = os.path.join(GP_RESOURCE_DIR, "effective_detection_rank_mapping.tsv")
    params:
        min_peptides = config["min_taxon_db_peptides"],
        resource_dir = GP_RESOURCE_DIR
    log:
        os.path.join(GP_RESOURCE_DIR, "logs/build_effective_detection_rank_db.log")
    container:
        config["containers"]["taxonkit"]
    shell:
        r"""
        set -euo pipefail
        LOG="{log}"
        MIN_PEP="{params.min_peptides}"
        LINEAGE_FILE="{input.taxid_family_map}"

        log_ts() {{ echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }}
        mkdir -p $(dirname "$LOG")

        log_ts "Building effective detection rank database from genome peptidotyping data"
        log_ts "Using pre-built lineage file: $LINEAGE_FILE"

        # ── Step 1: count peptides per family at each rank ───────────────────────
        COUNTS_FILE="{params.resource_dir}/family_rank_peptide_counts.tsv"
        echo -e "family_taxid\trank\tn_peptides" > "$COUNTS_FILE"

        # Family-level peptides: lca_il IS the family; count directly.
        tail -n +2 {input.family_tsv} | cut -f4 | awk '
            NF && $1!="" {{ count[$1]++ }}
            END {{ for (k in count) print k"\tfamily\t"count[k] }}
        ' >> "$COUNTS_FILE"

        # Genus-level peptides: map genus → family via lineage file.
        tail -n +2 {input.genus_tsv} | cut -f4 | awk -v lineage="$LINEAGE_FILE" '
            BEGIN {{ while((getline < lineage)>0) fam[$1]=$3 }}
            NF && $1!="" && $1 in fam {{ sum[fam[$1]]++ }}
            END {{ for (f in sum) print f"\tgenus\t"sum[f] }}
        ' >> "$COUNTS_FILE"

        # Species/strain-level peptides: map species → family via lineage file.
        tail -n +2 {input.species_tsv} | cut -f4 | awk -v lineage="$LINEAGE_FILE" '
            BEGIN {{ while((getline < lineage)>0) fam[$1]=$3 }}
            NF && $1!="" && $1 in fam {{ sum[fam[$1]]++ }}
            END {{ for (f in sum) print f"\tspecies_strain\t"sum[f] }}
        ' >> "$COUNTS_FILE"

        log_ts "Peptide counts per family/rank computed"

        # ── Step 2: determine effective detection rank per family ─────────────────
        MAPPING="{output.rank_mapping}"
        echo -e "family_taxid\teffective_rank\tn_peptides" > "$MAPPING"

        awk -F'\t' -v min="$MIN_PEP" '
            NR==1 {{ next }}
            {{
                fam=$1; rank=$2; n=$3
                count[fam"\t"rank] = n
                families[fam] = 1
            }}
            END {{
                for (fam in families) {{
                    fam_key = fam"\tfamily"
                    gen_key = fam"\tgenus"
                    sps_key = fam"\tspecies_strain"
                    if (fam_key in count && count[fam_key]+0 >= min+0) {{
                        print fam"\tfamily\t"count[fam_key]
                    }} else if (gen_key in count && count[gen_key]+0 >= min+0) {{
                        print fam"\tgenus\t"count[gen_key]
                    }} else if (sps_key in count && count[sps_key]+0 >= min+0) {{
                        print fam"\tspecies_strain\t"count[sps_key]
                    }}
                }}
            }}
        ' "$COUNTS_FILE" >> "$MAPPING"

        log_ts "Effective detection rank mapping: $(tail -n +2 "$MAPPING" | wc -l) families assigned"
        log_ts "  family rank:         $(awk -F'\t' '$2=="family"' "$MAPPING" | wc -l)"
        log_ts "  genus rank:          $(awk -F'\t' '$2=="genus"' "$MAPPING" | wc -l)"
        log_ts "  species_strain rank: $(awk -F'\t' '$2=="species_strain"' "$MAPPING" | wc -l)"

        # ── Step 3: build the first-pass FASTA ────────────────────────────────────
        > {output.first_pass_fasta}

        RANK_LOOKUP="{params.resource_dir}/effective_rank_lookup.tsv"
        tail -n +2 "$MAPPING" | awk -F'\t' '{{print $1"\t"$2}}' > "$RANK_LOOKUP"

        # Family-level peptides: include where family is assigned family rank
        awk -F'\t' '
            NR==FNR {{
                if ($2=="family") fam_families[$1]=1
                next
            }}
            FNR==1 {{ next }}
            $4 in fam_families {{
                fam_taxid = $4
                header = $10 " FAM=" fam_taxid
                print ">" header
                print $2
            }}
        ' "$RANK_LOOKUP" {input.family_tsv} >> {output.first_pass_fasta}

        # Genus-level peptides: include where parent family is assigned genus rank
        awk -F'\t' '
            ARGIND==1 {{ fam[$1]=$3; next }}
            ARGIND==2 {{ if($2=="genus") genus_fams[$1]=1; next }}
            FNR==1 {{ next }}
            {{
                lca_il=$4
                if (lca_il in fam && fam[lca_il] in genus_fams) {{
                    fam_taxid = fam[lca_il]
                    header = $10 " FAM=" fam_taxid
                    print ">" header
                    print $2
                }}
            }}
        ' "$LINEAGE_FILE" "$RANK_LOOKUP" {input.genus_tsv} >> {output.first_pass_fasta}

        # Species/strain-level peptides: include where parent family is assigned species_strain rank
        awk -F'\t' '
            ARGIND==1 {{ fam[$1]=$3; next }}
            ARGIND==2 {{ if($2=="species_strain") sps_fams[$1]=1; next }}
            FNR==1 {{ next }}
            {{
                lca_il=$4
                if (lca_il in fam && fam[lca_il] in sps_fams) {{
                    fam_taxid = fam[lca_il]
                    header = $10 " FAM=" fam_taxid
                    print ">" header
                    print $2
                }}
            }}
        ' "$LINEAGE_FILE" "$RANK_LOOKUP" {input.species_tsv} >> {output.first_pass_fasta}

        TOTAL=$(grep -c "^>" {output.first_pass_fasta} || true)
        log_ts "Effective first-pass database built: $TOTAL entries"
        """

################################################################################
# Phase D: Two-Pass DIA-NN Search
################################################################################

rule perform_genome_peptidotyping_first_pass_search:
    input:
        raw_files_dir = os.path.join(EXPERIMENT_DIR, "input/ms_files"),
        fasta = os.path.join(GP_RESOURCE_DIR, "effective_first_pass_database.fasta"),
        config_file = "config/peptidotyping_infinidia.cfg"
    output:
        first_pass_diann_parquet = os.path.join(GP_RESOURCE_DIR, "first_pass_diann.parquet"),
        first_pass_diann_protein_description = os.path.join(GP_RESOURCE_DIR, "first_pass_diann.protein_description.tsv")
    params:
        out_prefix = os.path.join(GP_RESOURCE_DIR, "first_pass_diann")
    log:
        os.path.join(RUN_DIR, "logs/genome_peptidotyping/first_pass_search.log")
    container:
        config["containers"]["diann"]
    threads: workflow.cores
    shell:
        """
        mkdir -p $(dirname {log})
        diann --cfg {input.config_file} \
            --fasta {input.fasta} \
            --out {params.out_prefix} \
            --dir {input.raw_files_dir} \
            --threads {threads} --verbose 1 >> {log} 2>&1
        """

rule infer_genome_peptidotyping_first_pass_presence:
    input:
        first_pass_diann = os.path.join(GP_RESOURCE_DIR, "first_pass_diann.parquet"),
        taxid_family_map = os.path.join(GP_RESOURCE_DIR, "taxid_to_family_genus.tsv")
    output:
        ncbi_taxonomy_id = os.path.join(GP_RESOURCE_DIR, "detected_family_taxa_ids.txt"),
        fdr_results      = os.path.join(GP_RESOURCE_DIR, "first_pass_fdr_results.tsv")
    log:
        os.path.join(RUN_DIR, "logs/genome_peptidotyping/infer_family_presence.log")
    container:
        config["containers"]["conduitr"]
    script:
        "../peptidotyping/scripts/infer_family_presence.R"

rule map_genome_peptidotyping_families_to_species:
    input:
        detected_families = os.path.join(GP_RESOURCE_DIR, "detected_family_taxa_ids.txt"),
        taxonomy = _taxonomy_input()
    output:
        families_to_species = os.path.join(GP_RESOURCE_DIR, "families_to_species.txt")
    log:
        os.path.join(RUN_DIR, "logs/genome_peptidotyping/map_families_to_species.log")
    container:
        config["containers"]["bakta"]
    script:
        "scripts/map_families_to_species.py"

rule generate_genome_peptidotyping_second_pass_db:
    input:
        species_tsv            = os.path.join(GP_RESOURCE_DIR, "species_strain_lca_filtered_peptides.tsv"),
        families_to_species    = os.path.join(GP_RESOURCE_DIR, "families_to_species.txt")
    output:
        second_pass_fasta = os.path.join(GP_RESOURCE_DIR, "second_pass_database.fasta")
    log:
        os.path.join(RUN_DIR, "logs/genome_peptidotyping/generate_second_pass_db.log")
    container:
        config["containers"]["taxonkit"]
    shell:
        r"""
        set -euo pipefail
        log_ts() {{ echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a {log}; }}
        mkdir -p $(dirname {log})

        log_ts "Building species allowlist from families_to_species.txt"

        ALLOWLIST=$(mktemp)
        awk 'NF && !/^[[:space:]]*$/' {input.families_to_species} \
            | tr -d '[:blank:]' | sort -u > "$ALLOWLIST"
        log_ts "Allowlist: $(wc -l < "$ALLOWLIST") species entries"

        log_ts "Filtering species_strain TSV to allowlist and writing FASTA"

        awk -F'\t' '
            ARGIND==1 {{ allowed[$1]=1; next }}
            FNR==1    {{ next }}
            $4 in allowed {{ print ">" $10; print $2 }}
        ' "$ALLOWLIST" {input.species_tsv} > {output.second_pass_fasta}

        rm -f "$ALLOWLIST"
        log_ts "Second-pass FASTA: $(grep -c "^>" {output.second_pass_fasta} || true) entries"
        """

rule perform_genome_peptidotyping_second_pass_search:
    input:
        raw_files_dir = os.path.join(EXPERIMENT_DIR, "input/ms_files"),
        fasta         = os.path.join(GP_RESOURCE_DIR, "second_pass_database.fasta"),
        config_file   = "config/peptidotyping_infinidia.cfg"
    output:
        second_pass_diann_parquet             = os.path.join(GP_RESOURCE_DIR, "second_pass_diann.parquet"),
        second_pass_diann_protein_description = os.path.join(GP_RESOURCE_DIR, "second_pass_diann.protein_description.tsv")
    params:
        out_prefix = os.path.join(GP_RESOURCE_DIR, "second_pass_diann")
    log:
        os.path.join(RUN_DIR, "logs/genome_peptidotyping/second_pass_search.log")
    container:
        config["containers"]["diann"]
    threads: workflow.cores
    shell:
        """
        mkdir -p $(dirname {log})
        diann --cfg {input.config_file} \
            --fasta {input.fasta} \
            --out {params.out_prefix} \
            --dir {input.raw_files_dir} \
            --threads {threads} --verbose 1 >> {log} 2>&1
        """

rule infer_genome_peptidotyping_second_pass_presence:
    input:
        second_pass_diann = os.path.join(GP_RESOURCE_DIR, "second_pass_diann.parquet")
    output:
        detected_species_strains = os.path.join(GP_RESOURCE_DIR, "detected_species_strain_taxa_ids.txt"),
        fdr_results              = os.path.join(GP_RESOURCE_DIR, "second_pass_fdr_results.tsv")
    log:
        os.path.join(RUN_DIR, "logs/genome_peptidotyping/infer_species_strain_presence.log")
    container:
        config["containers"]["conduitr"]
    script:
        "../peptidotyping/scripts/infer_species_strain_presence.R"

################################################################################
# Phase E: Genome Selection Checkpoint
################################################################################
# Maps detected species/strain taxa back to input genomes.
# This is a checkpoint — downstream MAGs rules re-evaluate which genomes to
# process based on this output.

checkpoint select_genomes_by_peptidotyping:
    input:
        detected_species = os.path.join(GP_RESOURCE_DIR, "detected_species_strain_taxa_ids.txt"),
        taxonomy = _taxonomy_input()
    output:
        detected_genomes = os.path.join(GP_RESOURCE_DIR, "detected_genomes.txt")
    log:
        os.path.join(RUN_DIR, "logs/genome_peptidotyping/select_genomes.log")
    container:
        config["containers"]["bakta"]
    script:
        "scripts/map_detected_species_to_genomes.py"
