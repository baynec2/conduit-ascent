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

EXPERIMENT_DIR = config["experiment_dir"]
RUN_DIR = config["run_dir"]
RAW_FILEPATHS = glob.glob(os.path.join(EXPERIMENT_DIR, "input/ms_files/*.raw"))

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
    container: config["containers"]["conduitr"]
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
rule perform_hapid_first_pass_search:
    input:
        raw_files_dir    = os.path.join(EXPERIMENT_DIR,"input/ms_files"),
        spectral_library = os.path.join(config["peptidotyping_resource_dir"],"hapid_peptidotyping.predicted.speclib"),
        fasta            = os.path.join(config["peptidotyping_resource_dir"],"hapid_peptidotyping_db.fasta"),
        config_file      = config["diann_library_search_base_config"]
    output:
        hapid_diann_parquet = os.path.join(RUN_DIR,"database_resources/unipept_hapid/hapid_first_pass_diann.parquet")
    log: os.path.join(RUN_DIR,"logs/unipept_hapid/perform_hapid_first_pass_search.log")
    container: config["containers"]["diann"]
    threads: workflow.cores
    shell:
        """
        diann --cfg {input.config_file} \
        --fasta {input.fasta} \
        --out  {RUN_DIR}/database_resources/unipept_hapid/hapid_first_pass_diann \
        --dir {input.raw_files_dir} \
        --lib {input.spectral_library} \
        --cut "" \
        --missed-cleavages 0 \
        --min-pep-len 5 \
        --max-pep-len 50 \
        --threads {threads} --verbose 1 >> {log} 2>&1
        """

################################################################################
# Inferring species/strain presence from HAPiID first-pass results
################################################################################
# Since hapid peptides already have species/strain-level LCA, presence is
# inferred directly at the species/strain level — no family→species mapping needed.
rule infer_species_presence:
    input:
        hapid_diann = os.path.join(RUN_DIR,"database_resources/unipept_hapid/hapid_first_pass_diann.parquet")
    output:
        ncbi_taxa_ids = os.path.join(RUN_DIR,"ncbi_taxa_ids.txt")
    params:
        threshold = config["presence_min_peptides"]
    log: os.path.join(RUN_DIR,"logs/unipept_hapid/infer_species_presence.log")
    container: config["containers"]["conduitr"]
    script: "scripts/infer_species_presence.R"
