################################################################################
# unipept_hapid search space module
################################################################################
# This module implements the HAPiID-inspired (Highly Abundant Protein
# Identification) search strategy using the Unipept LCA-assigned peptide
# database. Rather than using family-level taxonomic specificity (see the
# `peptidotyping` method), this approach filters for peptides whose LCA is at
# species/strain level AND whose functional annotation contains GO terms
# associated with highly abundant proteins (ribosomes, translation machinery).
#
# Because the peptides already have species/strain-level LCA, the pipeline is
# simpler than peptidotyping: there is no family→species mapping step. The
# first-pass search results are used directly to identify present species/strains
# which are then handed off as ncbi_taxa_ids.txt to the ncbi_taxonomy_id workflow.
#
# Note: build_sequence_index (which produces sequences.tsv.lz4 + taxons.tsv.lz4)
# is shared with the peptidotyping module and lives in
# modules/search_space/_shared/unipept_resources.smk.
#
# Reference: HAPiID — https://pmc.ncbi.nlm.nih.gov/articles/PMC8017886/
################################################################################

import glob
import os
import sys

EXPERIMENT_DIR = config["experiment_dir"]
RUN_DIR = config["run_dir"]

sys.path.insert(0, os.path.join(workflow.basedir, "modules", "_shared"))
from diann_staging import (
    list_raw_files,
    list_samples,
    raw_path_for_sample,
    stage3_symlink_commands,
)

RAW_FILEPATHS = list_raw_files(EXPERIMENT_DIR)
SAMPLES = list_samples(EXPERIMENT_DIR)

UH_OUT = os.path.join(RUN_DIR, "database_resources/unipept_hapid")
UH_QUANTS = os.path.join(UH_OUT, "first_pass_quant_files")

################################################################################
# Generating the HAPiID-style peptide database
################################################################################
# Peptides filtered to species/strain LCA AND annotated with GO terms for
# highly abundant proteins (ribosome, translation elongation factors).
# GO term selection follows the HAPiID rationale: if these proteins are absent
# from DIA-NN results, detection of lower-abundance proteins is implausible.
#
# Default GO terms:
#   GO:0005840  ribosome              ~22M annotations
#   GO:0006412  translation           ~23M annotations
#   GO:0003746  translation elongation factor activity  ~438K annotations
rule generate_hapid_database:
    input:
        sequences = os.path.join(config["peptidotyping_resource_dir"],"sequences.tsv.lz4"),
        taxons    = os.path.join(config["peptidotyping_resource_dir"],"taxons.tsv.lz4")
    output:
        lca_filtered_taxa = os.path.join(config["peptidotyping_resource_dir"],"hapid_lca_filtered_peptides.tsv"),
        hapid_fasta       = os.path.join(config["peptidotyping_resource_dir"],"hapid_peptidotyping_db.fasta")
    params:
        taxon_ranks_str = "species,strain",
        go_terms        = "GO:0005840,GO:0006412,GO:0003746"
    # Needs lz4 + awk + bash. Older conduitr:alpha had lz4; the rebuilt
    # conduitr (:f0dbc03 / :0ae9adc) dropped it. umgap has lz4 since it
    # produces these .lz4 indices, so use it here until conduitr regains lz4.
    container: config["containers"]["umgap"]
    log: os.path.join(config["peptidotyping_resource_dir"],"logs/generate_hapid_database.log")
    shell:
        r"""
        set -euo pipefail

        SEQUENCES_FILE="{input.sequences}"
        TAXONS_FILE="{input.taxons}"
        LCA_FILTERED_TAXA="{output.lca_filtered_taxa}"
        OUTPUT_FASTA="{output.hapid_fasta}"
        TAXON_RANKS="{params.taxon_ranks_str}"
        GO_TERMS="{params.go_terms}"
        LOG_FILE="{log}"

        log_with_timestamp() {{
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
        }}

        log_with_timestamp "Starting generate_hapid_database"
        log_with_timestamp "Taxon ranks: $TAXON_RANKS"
        log_with_timestamp "GO terms:    $GO_TERMS"

        RANK_PATTERN=$(echo "$TAXON_RANKS" | tr ',' '|')
        GO_PATTERN=$(echo "$GO_TERMS" | tr ',' '|')

        echo "$GO_TERMS" | tr ',' '\n' > /tmp/go_terms_list_$$.txt

        > "$LCA_FILTERED_TAXA"
        > "$OUTPUT_FASTA"
        echo -e "id\tsequence\tlca\tlca_il\tfa\tfa_il\tname\trank\tparent_id\tfasta_header" > "$LCA_FILTERED_TAXA"

        log_with_timestamp "Step 1: Counting eligible taxons"
        TAXON_COUNT=$(lz4 -d -c "$TAXONS_FILE" | awk -F'\t' -v rank_pattern="$RANK_PATTERN" '
            $3 ~ "^(" rank_pattern ")$" {{ count++ }}
            END {{ print count+0 }}
        ')
        log_with_timestamp "Taxons at target ranks: $TAXON_COUNT"

        log_with_timestamp "Step 2: Single-pass filter by taxonomy AND GO terms"
        awk -F'\t' \
            -v rank_pattern="$RANK_PATTERN" \
            -v go_pattern="$GO_PATTERN" \
            -v go_terms_file="/tmp/go_terms_list_$$.txt" \
            -v tsv_file="$LCA_FILTERED_TAXA" \
            -v fasta_file="$OUTPUT_FASTA" \
            -v log_file="$LOG_FILE" \
            '
            BEGIN {{
                go_term_count = 0
                while ((getline go_term < go_terms_file) > 0) {{
                    go_term_count++
                    go_terms_list[go_term_count] = go_term
                    go_term_peptide_count[go_term] = 0
                }}
                close(go_terms_file)
            }}
            FNR == NR {{
                if ($3 ~ "^(" rank_pattern ")$") {{
                    taxon_name[$1] = $2
                    taxon_rank[$1] = $3
                    taxon_parent[$1] = $4
                }}
                next
            }}
            {{
                total++
                seq_id = $1; sequence = $2; lca = $3; lca_il = $4; fa = $5; fa_il = $6

                taxonomy_match = (lca_il in taxon_name)
                go_match = 0
                if (fa_il != "" && fa_il != "\\N" && taxonomy_match) {{
                    for (i = 1; i <= go_term_count; i++) {{
                        go_term = go_terms_list[i]
                        if (match(fa_il, "\"" go_term "\"")) {{
                            go_match = 1
                            go_term_peptide_count[go_term]++
                            go_term_taxa[go_term SUBSEP lca_il] = 1
                        }}
                    }}
                }}

                if (taxonomy_match && go_match) {{
                    kept++
                    name = taxon_name[lca_il]; rank = taxon_rank[lca_il]; parent_id = taxon_parent[lca_il]
                    name_for_header = name; gsub(/ /, "-", name_for_header)
                    fasta_header = "umgap|" seq_id "|" lca_il " " rank "_" name_for_header " OS=" name " OX=" lca_il " RK=" rank " PT=" parent_id
                    print seq_id "\t" sequence "\t" lca "\t" lca_il "\t" fa "\t" fa_il "\t" name "\t" rank "\t" parent_id "\t" fasta_header >> tsv_file
                    print ">" fasta_header >> fasta_file
                    print sequence >> fasta_file
                }}

                if (total % 10000000 == 0) {{
                    printf "[%s] Processed %d sequences, kept %d (%.2f%%)\n", strftime("%Y-%m-%d %H:%M:%S"), total, kept, (kept/total)*100 >> log_file
                    close(log_file)
                }}
            }}
            END {{
                printf "[%s] Filtered to %d peptides (from %d total, %.2f%%)\n", strftime("%Y-%m-%d %H:%M:%S"), kept, total, (kept/total)*100 >> log_file
                printf "\n[%s] GO Term Statistics:\nGO_Term\tPeptides\tUnique_Taxa\n", strftime("%Y-%m-%d %H:%M:%S") >> log_file
                for (i = 1; i <= go_term_count; i++) {{
                    go_term = go_terms_list[i]
                    taxa_count = 0
                    for (key in go_term_taxa) {{
                        split(key, arr, SUBSEP)
                        if (arr[1] == go_term) taxa_count++
                    }}
                    printf "%s\t%d\t%d\n", go_term, go_term_peptide_count[go_term], taxa_count >> log_file
                }}
                close(log_file)
            }}
            ' <(lz4 -d -c "$TAXONS_FILE") <(lz4 -d -c "$SEQUENCES_FILE")

        rm -f /tmp/go_terms_list_$$.txt
        log_with_timestamp "Completed generate_hapid_database"
        """

################################################################################
# Generating the HAPiID Spectral Library
################################################################################
rule generate_hapid_spectral_library:
    input:
        fasta       = os.path.join(config["peptidotyping_resource_dir"],"hapid_peptidotyping_db.fasta"),
        config_file = config["diann_spectral_library_base_config"]
    output:
        os.path.join(config["peptidotyping_resource_dir"],"hapid_peptidotyping.predicted.speclib")
    # DIA-NN appends .predicted.speclib to the --out-lib path
    container: config["containers"]["diann"]
    log: os.path.join(config["peptidotyping_resource_dir"],"logs/generate_hapid_spectral_library.log")
    threads: workflow.cores
    shell:
        """
        diann --cfg {input.config_file} \
        --fasta {input.fasta} \
        --threads {threads} \
        --out-lib {config[peptidotyping_resource_dir]}hapid_peptidotyping \
        --cut "" \
        --missed-cleavages 0 \
        --min-pep-len 5 \
        --max-pep-len 50 \
        --species-ids >> {log} 2>&1
        """

################################################################################
# Performing the HAPiID First Pass Search
################################################################################
# Standard: 3-stage split. InfinDIA: monolithic. See modules/diann/diann.smk for rationale.
if config.get("unipept_hapid_search_mode", "standard") == "standard":

    rule unipept_hapid_build_empirical_lib:
        input:
            raw_files_dir    = os.path.join(EXPERIMENT_DIR,"input/ms_files"),
            spectral_library = os.path.join(config["peptidotyping_resource_dir"],"hapid_peptidotyping.predicted.speclib"),
            fasta            = os.path.join(config["peptidotyping_resource_dir"],"hapid_peptidotyping_db.fasta"),
            config_file      = os.path.join(RUN_DIR,"config/diann_library_search_base.cfg")
        output:
            empirical_lib = os.path.join(UH_OUT, "hapid_first_pass_empirical.parquet")
        params:
            out_lib = os.path.join(UH_OUT, "hapid_first_pass_empirical"),
            tmpdir = os.path.join(UH_OUT, "build_empirical_quant_files"),
        log: os.path.join(RUN_DIR,"logs/unipept_hapid/build_empirical_lib.log")
        container: config["containers"]["diann"]
        threads: workflow.cores
        shell:
            """
            mkdir -p $(dirname {log}) $(dirname {output.empirical_lib})
            rm -rf {params.tmpdir} && mkdir -p {params.tmpdir}
            diann --cfg {input.config_file} \
            --fasta {input.fasta} \
            --dir {input.raw_files_dir} \
            --temp {params.tmpdir} \
            --lib {input.spectral_library} \
            --gen-spec-lib \
            --rt-profiling \
            --out-lib {params.out_lib} \
            --out {params.tmpdir}/report \
            --cut "" \
            --missed-cleavages 0 \
            --min-pep-len 5 \
            --max-pep-len 50 \
            --threads {threads} --verbose 1 >> {log} 2>&1
            """

    rule unipept_hapid_search_one_raw:
        input:
            empirical_lib = os.path.join(UH_OUT, "hapid_first_pass_empirical.parquet"),
            fasta = os.path.join(config["peptidotyping_resource_dir"],"hapid_peptidotyping_db.fasta"),
            config_file = os.path.join(RUN_DIR,"config/diann_library_search_base.cfg"),
            raw = lambda w: raw_path_for_sample(EXPERIMENT_DIR, w.sample)
        output:
            quant = os.path.join(UH_QUANTS, "{sample}.quant")
        params:
            tmpdir = lambda w: os.path.join(UH_OUT, "first_pass_quant_tmp", w.sample)
        log: os.path.join(RUN_DIR,"logs/unipept_hapid/search_one_raw.{sample}.log")
        container: config["containers"]["diann"]
        threads: min(8, workflow.cores)
        shell:
            """
            mkdir -p $(dirname {log})
            rm -rf {params.tmpdir} && mkdir -p {params.tmpdir} $(dirname {output.quant})
            diann --cfg {input.config_file} \
            --f {input.raw} \
            --lib {input.empirical_lib} \
            --fasta {input.fasta} \
            --temp {params.tmpdir} \
            --out {params.tmpdir}/per_run_report \
            --cut "" \
            --missed-cleavages 0 \
            --min-pep-len 5 \
            --max-pep-len 50 \
            --threads {threads} --verbose 1 >> {log} 2>&1
            mv {params.tmpdir}/*.quant {output.quant}
            rm -rf {params.tmpdir}
            """

    rule unipept_hapid_combine:
        input:
            quants = expand(
                os.path.join(UH_QUANTS, "{sample}.quant"),
                sample=SAMPLES
            ),
            empirical_lib = os.path.join(UH_OUT, "hapid_first_pass_empirical.parquet"),
            fasta = os.path.join(config["peptidotyping_resource_dir"],"hapid_peptidotyping_db.fasta"),
            config_file = os.path.join(RUN_DIR,"config/diann_library_search_base.cfg"),
            raw_files_dir = os.path.join(EXPERIMENT_DIR,"input/ms_files")
        output:
            hapid_diann_parquet = os.path.join(UH_OUT, "hapid_first_pass_diann.parquet")
        params:
            tmpdir = os.path.join(UH_OUT, "first_pass_combine_tmp"),
            out_prefix = os.path.join(UH_OUT, "hapid_first_pass_diann"),
            symlink_cmds = stage3_symlink_commands(
                os.path.join(UH_OUT, "first_pass_combine_tmp"),
                RAW_FILEPATHS,
                UH_QUANTS
            )
        log: os.path.join(RUN_DIR,"logs/unipept_hapid/combine.log")
        container: config["containers"]["diann"]
        threads: workflow.cores
        shell:
            """
            mkdir -p $(dirname {log})
            rm -rf {params.tmpdir} && mkdir -p {params.tmpdir}
            {params.symlink_cmds}
            diann --cfg {input.config_file} \
            --dir {input.raw_files_dir} \
            --lib {input.empirical_lib} \
            --fasta {input.fasta} \
            --temp {params.tmpdir} \
            --use-quant \
            --out {params.out_prefix} \
            --cut "" \
            --missed-cleavages 0 \
            --min-pep-len 5 \
            --max-pep-len 50 \
            --threads {threads} --verbose 1 >> {log} 2>&1
            rm -rf {params.tmpdir}
            """

else:  # infinidia — monolithic

    rule unipept_hapid_monolithic:
        input:
            raw_files_dir = os.path.join(EXPERIMENT_DIR,"input/ms_files"),
            fasta         = os.path.join(config["peptidotyping_resource_dir"],"hapid_peptidotyping_db.fasta"),
            config_file   = os.path.join(RUN_DIR,"config/hapid_infinidia.cfg")
        output:
            hapid_diann_parquet = os.path.join(UH_OUT, "hapid_first_pass_diann.parquet")
        params:
            out_prefix = os.path.join(UH_OUT, "hapid_first_pass_diann"),
            tmpdir = os.path.join(UH_OUT, "monolithic_quant_files"),
        log: os.path.join(RUN_DIR,"logs/unipept_hapid/monolithic.log")
        container: config["containers"]["diann"]
        threads: workflow.cores
        shell:
            """
            mkdir -p $(dirname {log}) $(dirname {output.hapid_diann_parquet})
            rm -rf {params.tmpdir} && mkdir -p {params.tmpdir}
            diann --cfg {input.config_file} \
            --fasta {input.fasta} \
            --dir {input.raw_files_dir} \
            --temp {params.tmpdir} \
            --pre-search --pre-filter \
            --gen-spec-lib \
            --rt-profiling \
            --out {params.out_prefix} \
            --cut "" \
            --missed-cleavages 0 \
            --min-pep-len 5 \
            --max-pep-len 50 \
            --threads {threads} --verbose 1 >> {log} 2>&1
            """

################################################################################
# Greedy HAPiID-style species/strain selection from first-pass results
################################################################################
# Mirrors the genome-based hapid module (modules/search_space/hapid):
#   1. Build a taxon → {spectrum_ids} dict (spectrum = Run||Precursor.Id).
#   2. Greedy set-cover: at each step pick the taxon covering the most
#      still-uncovered spectra; emit a TSV ranked by cumulative coverage.
#   3. Take the smallest prefix whose cumulative_pct ≥ hapid_percent_spectra
#      and write it to ncbi_taxa_ids.txt for the ncbi_taxonomy_id workflow.
#
# The greedy script lives in modules/search_space/_shared/scripts/ and is
# reused as-is across the genome-based hapid module and this one.
rule build_taxon_spectrum_mapping:
    input:
        parquet = os.path.join(RUN_DIR, "database_resources/unipept_hapid/hapid_first_pass_diann.parquet")
    output:
        os.path.join(RUN_DIR, "database_resources/unipept_hapid/taxon2spectrum_dic.json")
    log: os.path.join(RUN_DIR, "logs/unipept_hapid/build_taxon_spectrum_mapping.log")
    container: config["containers"]["fraggenescan_hmmer"]
    script: "scripts/build_taxon_spectrum_mapping.py"


checkpoint run_greedy_taxon_selection:
    input:
        os.path.join(RUN_DIR, "database_resources/unipept_hapid/taxon2spectrum_dic.json")
    output:
        os.path.join(RUN_DIR, "database_resources/unipept_hapid/unipept_hapid_greedy_selection_raw.tsv")
    log: os.path.join(RUN_DIR, "logs/unipept_hapid/greedy_taxon_selection.log")
    container: config["containers"]["fraggenescan_hmmer"]
    shell:
        """
        mkdir -p $(dirname {log})
        python {workflow.basedir}/modules/search_space/_shared/scripts/coverAllSpectra_greedy.py \
            {input} {output} > {log} 2>&1
        """


# Enrich the greedy TSV with NCBI scientific names. The `genome` column is an
# NCBI taxid (lca_il); the name is already in hapid_lca_filtered_peptides.tsv
# from the database-build step, so we just dedupe + join — no taxonkit needed.
rule annotate_unipept_hapid_greedy_selection:
    input:
        raw     = os.path.join(RUN_DIR, "database_resources/unipept_hapid/unipept_hapid_greedy_selection_raw.tsv"),
        lca_tsv = os.path.join(config["peptidotyping_resource_dir"], "hapid_lca_filtered_peptides.tsv")
    output:
        os.path.join(RUN_DIR, "database_resources/unipept_hapid/unipept_hapid_greedy_selection.tsv")
    log: os.path.join(RUN_DIR, "logs/unipept_hapid/annotate_greedy_selection.log")
    run:
        import pandas as pd
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        raw = pd.read_csv(input.raw, sep="\t", dtype={"genome": str})
        names = (
            pd.read_csv(input.lca_tsv, sep="\t", usecols=["lca_il", "name"], dtype={"lca_il": str})
              .drop_duplicates(subset=["lca_il"])
              .rename(columns={"lca_il": "genome", "name": "taxon_name"})
        )
        merged = raw.merge(names, on="genome", how="left")
        merged = merged[["genome", "taxon_name", "nSpectraCovered", "cumulative_pct"]]
        merged.to_csv(output[0], sep="\t", index=False)
        with open(log[0], "w") as lh:
            n_missing = merged["taxon_name"].isna().sum()
            lh.write(
                f"Annotated {len(merged)} taxa with scientific names "
                f"({n_missing} missing).\n"
            )


# Apply the hapid_percent_spectra cutoff and emit ncbi_taxa_ids.txt.
# The shared greedy script writes its key column literally as `genome`; here
# the keys are NCBI taxon IDs, and we rename only at this output boundary.
rule unipept_hapid_filter_selected_taxa:
    input:
        os.path.join(RUN_DIR, "database_resources/unipept_hapid/unipept_hapid_greedy_selection.tsv")
    output:
        os.path.join(RUN_DIR, "ncbi_taxa_ids.txt")
    params:
        pct = config.get("hapid_percent_spectra", 80)
    log: os.path.join(RUN_DIR, "logs/unipept_hapid/filter_selected_taxa.log")
    run:
        import pandas as pd
        df = pd.read_csv(input[0], sep="\t")
        above = df[df["cumulative_pct"] >= params.pct]
        cutoff = (above.index[0] + 1) if not above.empty else len(df)
        selected = df["genome"].tolist()[:cutoff]
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        with open(output[0], "w") as fh:
            fh.write("ncbi_taxonomy_id\n")
            for t in selected:
                fh.write(f"{t}\n")
        with open(log[0], "w") as lh:
            lh.write(
                f"Selected {len(selected)} of {len(df)} taxa "
                f"(cumulative_pct ≥ {params.pct}); written to {output[0]}\n"
            )
