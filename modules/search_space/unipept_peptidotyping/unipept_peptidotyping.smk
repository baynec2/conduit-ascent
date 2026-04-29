import glob
import os
EXPERIMENT_DIR = config["experiment_dir"]
RUN_DIR = config["run_dir"]
RAW_FILEPATHS = glob.glob(os.path.join(EXPERIMENT_DIR, "input/ms_files/*.raw"))

################################################################################
# Determining Version of Uniprotkb that is being used in experiment
################################################################################
# Note: build_sequence_index (which produces sequences.tsv.lz4 + taxons.tsv.lz4)
# is shared with the unipept_hapid module and lives in
# modules/search_space/_shared/unipept_resources.smk.
rule check_sequence_index_version:
    input:
        relnotes = os.path.join(config["peptidotyping_resource_dir"],"relnotes.txt")
    output:
        version_file = os.path.join(config["peptidotyping_resource_dir"],".version")
    params:
        outdir = config["peptidotyping_resource_dir"]
    log:
        os.path.join(RUN_DIR,"logs/search_space/check_sequence_index_version.log")
    shell:
        """
        # Extract version from local relnotes.txt
        LOCAL_VERSION=$(head -1 {input} | grep -o 'Release [0-9_]*' | cut -d' ' -f2)
        
        # Get current version from UniProt
        CURRENT_VERSION=$(curl -s https://ftp.uniprot.org/pub/databases/uniprot/relnotes.txt | head -1 | grep -o 'Release [0-9_]*' | cut -d' ' -f2)
        
        # Write results to log
        echo "Used version: $LOCAL_VERSION" > {log}
        echo "Current version: $CURRENT_VERSION" >> {log}
        
        if [ "$LOCAL_VERSION" = "$CURRENT_VERSION" ]; then
            echo "Sequence index is up to date." >> {log}
        else
            echo "Sequence index is NOT up to date." >> {log}
        fi
        """
################################################################################
# Extracting Metrics from Peptidotyping Resources
################################################################################
# This rule extracts metrics from the peptidotyping resources to help understand
# the composition of the database.
rule extract_peptidotyping_resource_metrics:
    input:
        sequences = os.path.join(config["peptidotyping_resource_dir"],"sequences.tsv.lz4"),
        taxons = os.path.join(config["peptidotyping_resource_dir"],"taxons.tsv.lz4")
    output:
        peptidotyping_resource_metrics = os.path.join(config["peptidotyping_resource_dir"],"peptidotyping_resource_metrics.tsv")
    container: config["containers"]["conduitr"]
    log:
        os.path.join(config["peptidotyping_resource_dir"],"logs/extract_peptidotyping_resource_metrics.log")
    script:
        "scripts/extract_peptidotyping_resource_metrics.R"

################################################################################
# Specifying the peptidotyping fasta databases to generate for the next rule.
################################################################################
# Map rank name (used in filenames) to comma-separated taxon ranks for generate_peptidotyping_db
PEPTIDOTYPING_RANK_CONFIG = {
    "species_strain": "species,strain",
    "genus": "genus",
    "family": "family",
}


################################################################################
# Generating peptidotyping databases
################################################################################
# From the file containing the peptides and their LCAS, we can generate a fasta file
# only containing proteotypic peptides for a given taxonomy. We can subsequently use 
# this as a first pass database to identify taxa that are likely present in the experiment.
# The idea here is to first use broad taxonomic ranks to identify taxa that are likely present in the experiment.
# In the process, many taxa should be excluded, which will reduce the number of their children taxa, a form of hierarchical filtering.
# This will allow us to filter the species level peptides to a more manageable number.
# I am currently not sure what level of taxonomic rank would perform the best for the first pass search, so we will have to test a few things. 
# I am inclined to think phylum would perform the best, as it is the broadest rank and will likely exclude the most taxa.
# But there are arguments to be made for class, order, family, etc.
# We will test phylum, class, and family.

rule generate_peptidotyping_db:
    input:
        sequences = os.path.join(config["peptidotyping_resource_dir"],"sequences.tsv.lz4"),
        taxons = os.path.join(config["peptidotyping_resource_dir"],"taxons.tsv.lz4")
    output:
        lca_filtered_taxa = os.path.join(config["peptidotyping_resource_dir"],"{rank}_lca_filtered_peptides.tsv"),
        first_pass_fasta = os.path.join(config["peptidotyping_resource_dir"],"{rank}_peptidotyping_db.fasta")
    params:
        taxon_ranks_str = lambda wildcards: PEPTIDOTYPING_RANK_CONFIG[wildcards.rank]
    container: config["containers"]["conduitr"]
    log: os.path.join(config["peptidotyping_resource_dir"],"logs/generate_peptidotyping_db_{rank}.log")
    shell:
        r"""
        set -euo pipefail
        
        SEQUENCES_FILE="{input.sequences}"
        TAXONS_FILE="{input.taxons}"
        LCA_FILTERED_TAXA="{output.lca_filtered_taxa}"
        OUTPUT_FASTA="{output.first_pass_fasta}"
        TAXON_RANKS="{params.taxon_ranks_str}"
        LOG_FILE="{log}"
        
        # Function to log with timestamp
        log_with_timestamp() {{
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
        }}
        
        log_with_timestamp "Starting generate_peptidotyping_db_awk"
        log_with_timestamp "Input sequences: $SEQUENCES_FILE"
        log_with_timestamp "Input taxons: $TAXONS_FILE"
        log_with_timestamp "Output TSV: $LCA_FILTERED_TAXA"
        log_with_timestamp "Output FASTA: $OUTPUT_FASTA"
        log_with_timestamp "Taxon ranks: $TAXON_RANKS"
        
        # Convert comma-separated ranks to awk pattern
        RANK_PATTERN=$(echo "$TAXON_RANKS" | tr ',' '|')
        
        # Create empty output files
        > "$LCA_FILTERED_TAXA"
        > "$OUTPUT_FASTA"
        
        # Write TSV header
        echo -e "id\tsequence\tlca\tlca_il\tfa\tfa_il\tname\trank\tparent_id\tfasta_header" > "$LCA_FILTERED_TAXA"
        
        log_with_timestamp "Step 1: Building taxonomy lookup (tiny, once)"
        
        # Count taxons first (for logging)
        TAXON_COUNT=$(lz4 -d -c "$TAXONS_FILE" | awk -F'\t' -v rank_pattern="$RANK_PATTERN" '
            $3 ~ "^(" rank_pattern ")$" {{
                count++
            }}
            END {{
                print count+0
            }}
        ')
        
        log_with_timestamp "Will load $TAXON_COUNT taxons at ranks: $TAXON_RANKS"
        
        log_with_timestamp "Step 2: Streaming sequences, filtering, and writing output (single pass)"
        
        # Single pass - stream sequences, filter, and write TSV + FASTA
        awk -F'\t' \
            -v rank_pattern="$RANK_PATTERN" \
            -v tsv_file="$LCA_FILTERED_TAXA" \
            -v fasta_file="$OUTPUT_FASTA" \
            -v log_file="$LOG_FILE" \
            '
            # Process first file (taxons): build lookup
            FNR == NR {{
                if ($3 ~ "^(" rank_pattern ")$") {{
                    taxon_name[$1] = $2
                    taxon_rank[$1] = $3
                    taxon_parent[$1] = $4
                }}
                next
            }}
            
            # Process second file (sequences): filter and write
            {{
                total++
                
                seq_id = $1
                sequence = $2
                lca = $3
                lca_il = $4
                fa = $5
                fa_il = $6
                
                if (lca_il in taxon_name) {{
                    kept++
                    
                    name = taxon_name[lca_il]
                    rank = taxon_rank[lca_il]
                    parent_id = taxon_parent[lca_il]
                    
                    name_for_header = name
                    gsub(/ /, "-", name_for_header)
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
                printf "[%s] Filtered to %d peptides with LCA at specified ranks (from %d total, %.2f%%)\n", strftime("%Y-%m-%d %H:%M:%S"), kept, total, (kept/total)*100 >> log_file
                close(log_file)
            }}
            ' <(lz4 -d -c "$TAXONS_FILE") <(lz4 -d -c "$SEQUENCES_FILE")
        
        log_with_timestamp "Completed generate_peptidotyping_db_awk"
        log_with_timestamp "TSV file: $LCA_FILTERED_TAXA"
        log_with_timestamp "FASTA file: $OUTPUT_FASTA"
        """
################################################################################
# Building the Effective Detection Rank Database
################################################################################
# For each family in the taxonomy, the effective detection rank is the finest
# LCA rank (family → genus → species/strain) at which ≥ min_taxon_db_peptides
# unique proteotypic peptides exist. This approach replaces the previous
# genus-fallback mechanism (which was broken because genus-level peptides were
# never generated) with a principled, multi-rank strategy.
#
# Algorithm:
#   1. Build a taxid → family_taxid lineage lookup via taxonkit.
#   2. For each family: aggregate peptide counts across ALL child taxa at each rank
#      (family, genus, species/strain), producing one total per (family, rank).
#   3. Assign effective rank: finest rank where the family's aggregate total
#      >= min_taxon_db_peptides (family → genus → species/strain).
#   4. Build the first-pass FASTA including ALL child taxa at the effective rank
#      (not just a single representative), so the family can be detected even when
#      no single child taxon has enough peptides on its own.
#   5. Add FAM=<family_taxid> to every FASTA header so infer_family_presence.R
#      can always recover the family from any hit regardless of the rank.
#   6. Output effective_detection_rank_mapping.tsv for use in inference.
rule build_effective_detection_rank_db:
    input:
        family_tsv  = os.path.join(config["peptidotyping_resource_dir"],"family_lca_filtered_peptides.tsv"),
        genus_tsv   = os.path.join(config["peptidotyping_resource_dir"],"genus_lca_filtered_peptides.tsv"),
        species_tsv = os.path.join(config["peptidotyping_resource_dir"],"species_strain_lca_filtered_peptides.tsv")
    output:
        first_pass_fasta    = os.path.join(config["peptidotyping_resource_dir"],"effective_first_pass_database.fasta"),
        rank_mapping        = os.path.join(config["peptidotyping_resource_dir"],"effective_detection_rank_mapping.tsv"),
        taxid_family_map    = os.path.join(config["peptidotyping_resource_dir"],"taxid_to_family_genus.tsv")
    params:
        min_peptides = config["min_taxon_db_peptides"]
    threads: workflow.cores
    container: config["containers"]["taxonkit"]
    log: os.path.join(config["peptidotyping_resource_dir"],"logs/build_effective_detection_rank_db.log")
    shell:
        r"""
        set -euo pipefail
        LOG="{log}"
        MIN_PEP="{params.min_peptides}"

        log_ts() {{ echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }}

        export TAXONKIT_DB="{config[taxonkit_db_dir]}"

        log_ts "Building taxonomy lineage lookup (family/genus/species hierarchy)"

        # ── Step 1: build taxid → parent_family mapping using taxonkit ──────────
        # We need to know: for each genus taxid, which family does it belong to?
        # For each species/strain taxid, which family and genus does it belong to?
        LINEAGE_FILE="{config[peptidotyping_resource_dir]}/taxid_to_family_genus.tsv"

        # Collect all unique lca_il values across all three databases
        {{
          tail -n +2 {input.family_tsv}  | cut -f4
          tail -n +2 {input.genus_tsv}   | cut -f4
          tail -n +2 {input.species_tsv} | cut -f4
        }} | awk 'NF && $1!="" && !seen[$1]++' \
        | taxonkit lineage -j {threads} -t -R -r \
        | awk -F'\t' '{{
            n=split($5,rk,";"); split($3,lt,";");
            fam=""; gen="";
            for(i=1;i<=n;i++) {{
                if(rk[i]=="family") fam=lt[i];
                if(rk[i]=="genus")  gen=lt[i];
            }}
            if (fam!="") print $1"\t"$4"\t"fam"\t"gen
        }}' > "$LINEAGE_FILE"

        log_ts "Lineage lookup built: $(wc -l < "$LINEAGE_FILE") entries"

        # ── Step 2: count peptides per family at each rank ───────────────────────
        # family_counts[family_taxid][rank] = n_peptides

        COUNTS_FILE="{config[peptidotyping_resource_dir]}/family_rank_peptide_counts.tsv"
        # One row per (family, rank): family_taxid | rank | total_n_peptides
        echo -e "family_taxid\trank\tn_peptides" > "$COUNTS_FILE"

        # Family-level peptides: lca_il IS the family taxid; count directly in awk.
        tail -n +2 {input.family_tsv} | cut -f4 | awk '
            NF && $1!="" {{ count[$1]++ }}
            END {{ for (k in count) print k"\tfamily\t"count[k] }}
        ' >> "$COUNTS_FILE"

        # Genus-level peptides: stream taxids and accumulate into parent family in one pass.
        tail -n +2 {input.genus_tsv} | cut -f4 | awk -v lineage="$LINEAGE_FILE" '
            BEGIN {{ while((getline < lineage)>0) fam[$1]=$3 }}
            NF && $1!="" && $1 in fam {{ sum[fam[$1]]++ }}
            END {{ for (f in sum) print f"\tgenus\t"sum[f] }}
        ' >> "$COUNTS_FILE"

        # Species/strain-level peptides: same single-pass approach.
        tail -n +2 {input.species_tsv} | cut -f4 | awk -v lineage="$LINEAGE_FILE" '
            BEGIN {{ while((getline < lineage)>0) fam[$1]=$3 }}
            NF && $1!="" && $1 in fam {{ sum[fam[$1]]++ }}
            END {{ for (f in sum) print f"\tspecies_strain\t"sum[f] }}
        ' >> "$COUNTS_FILE"

        log_ts "Peptide counts per family/rank computed"

        # ── Step 3: determine effective detection rank per family ─────────────────
        # Priority: family > genus > species_strain (finest available with >= min_peptides)

        MAPPING="{output.rank_mapping}"
        # Columns: family_taxid | effective_rank | n_peptides
        echo -e "family_taxid\teffective_rank\tn_peptides" > "$MAPPING"

        # COUNTS_FILE has exactly one row per (family, rank) with an aggregate total,
        # so we just look up each family's count at each rank directly. Pick the finest
        # rank whose aggregate peptide total meets the threshold (family → genus → species/strain).
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
                    # else: family not reliably detectable; excluded from first-pass
                }}
            }}
        ' "$COUNTS_FILE" >> "$MAPPING"

        log_ts "Effective detection rank mapping: $(tail -n +2 "$MAPPING" | wc -l) families assigned"
        log_ts "  family rank:         $(awk -F'\t' '$2=="family"' "$MAPPING" | wc -l)"
        log_ts "  genus rank:          $(awk -F'\t' '$2=="genus"' "$MAPPING" | wc -l)"
        log_ts "  species_strain rank: $(awk -F'\t' '$2=="species_strain"' "$MAPPING" | wc -l)"

        # ── Step 4: build the first-pass FASTA ────────────────────────────────────
        # For each family in the mapping, include entries from the appropriate TSV
        # and add FAM=<family_taxid> to every header.

        > {output.first_pass_fasta}

        # Build a lookup of effective rank for each family
        RANK_LOOKUP="{config[peptidotyping_resource_dir]}/effective_rank_lookup.tsv"
        tail -n +2 "$MAPPING" | awk -F'\t' '{{print $1"\t"$2}}' > "$RANK_LOOKUP"

        # Process family-level peptides: include entries where family is in rank_lookup with rank==family
        # Add FAM= tag to the header (lca_il IS the family taxid for family-level entries)
        awk -F'\t' '
            NR==FNR {{
                if ($2=="family") fam_families[$1]=1
                next
            }}
            FNR==1 {{ next }}
            $4 in fam_families {{
                fam_taxid = $4
                # Rewrite fasta_header (col 10) to add FAM= tag
                header = $10 " FAM=" fam_taxid
                print ">" header
                print $2
            }}
        ' "$RANK_LOOKUP" {input.family_tsv} >> {output.first_pass_fasta}

        # Process genus-level peptides: include entries where parent family is in rank_lookup with rank==genus
        # Build family→genus lookup from lineage file
        awk -F'\t' '
            NR==FNR {{
                if ($2=="genus") genus_families[$1]=1
                next
            }}
            FNR==1 {{ next }}
            {{
                # lca_il ($4) is a genus taxid; look up its family from lineage
                # We need lineage info to find the family_taxid for this genus entry
                # lineage_file has: taxid, rank, family_taxid, genus_taxid
                # For genus-level entries in genus_tsv, lca_il IS the genus; lineage has fam[$4]
                # So we need to join genus_tsv with lineage on lca_il
            }}
        ' "$RANK_LOOKUP" {input.genus_tsv} > /dev/null  # placeholder: actual join below

        # More efficient: single awk pass joining genus_tsv with lineage + rank_lookup
        awk -F'\t' '
            # File 1: lineage (taxid, rank, family_taxid, genus_taxid)
            ARGIND==1 {{ fam[$1]=$3; next }}
            # File 2: rank_lookup (family_taxid, effective_rank)
            ARGIND==2 {{ if($2=="genus") genus_fams[$1]=1; next }}
            # File 3: genus_tsv (id, sequence, lca, lca_il, fa, fa_il, name, rank, parent_id, fasta_header)
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

        # Process species/strain-level peptides: include entries for families assigned species_strain rank
        awk -F'\t' '
            # File 1: lineage
            ARGIND==1 {{ fam[$1]=$3; next }}
            # File 2: rank_lookup
            ARGIND==2 {{ if($2=="species_strain") sps_fams[$1]=1; next }}
            # File 3: species_tsv
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
# Performing the First Pass Search, Find Broad Taxonomic Levels
################################################################################
# Searching our first pass spectral library with Diann on infinidia mode. 
# The database searched is mostly family level peptides, but in the case that no
# family level peptides for a family exists, fallback to genus or species/strain
# specific if they do not exist (as defined in uild_effective_detection_rank_db).
rule perform_first_pass_search:
    input:
        raw_files_dir = os.path.join(EXPERIMENT_DIR,"input/ms_files"),
        fasta = os.path.join(config["peptidotyping_resource_dir"],"effective_first_pass_database.fasta"),
        config_file = "config/peptidotyping_infinidia.cfg"
    output:
        first_pass_diann_parquet = os.path.join(RUN_DIR,"database_resources/peptidotyping/first_pass_diann.parquet"),
        first_pass_diann_protein_description =  os.path.join(RUN_DIR,"database_resources/peptidotyping/first_pass_diann.protein_description.tsv")
    log: os.path.join(RUN_DIR,"logs/peptidotyping/perfrom_first_pass_search.log")
    container:
        config["containers"]["diann"]
    threads: workflow.cores
    shell:
        """
        diann --cfg {input.config_file} \
        --fasta {input.fasta} \
        --out  {RUN_DIR}/database_resources/peptidotyping/first_pass_diann \
        --dir {input.raw_files_dir} \
        --threads {threads} --verbose 1 >> {log} 2>&1
        """
################################################################################
# Determine what families are present based on the first pass results
################################################################################
# Here we will use a taxonomic target decoy strategy to determine which taxonomic
# levels are likely present. 
rule infer_first_pass_presence:
    input:
        first_pass_diann  = os.path.join(RUN_DIR,"database_resources/peptidotyping/first_pass_diann.parquet"),
        taxid_family_map  = os.path.join(config["peptidotyping_resource_dir"],"taxid_to_family_genus.tsv")
    output:
        ncbi_taxonomy_id = os.path.join(RUN_DIR,"database_resources/peptidotyping/detected_family_taxa_ids.txt"),
        fdr_results      = os.path.join(RUN_DIR,"database_resources/peptidotyping/first_pass_fdr_results.tsv")
    log: os.path.join(RUN_DIR,"logs/peptidotyping/infer_family_presence.log")
    container: config["containers"]["conduitr"]
    script: "scripts/infer_family_presence.R"

################################################################################
# Determine which species and strains are descendants of the taxa detected in 
# the first pass search.
################################################################################
# We will use the species/strain specific peptides for a second pass search to 
# determine what species/strains are present. 
rule map_first_pass_detected_taxa_to_species_strains:
    input:
        ncbi_taxonomy_ids = os.path.join(RUN_DIR,"database_resources/peptidotyping/detected_family_taxa_ids.txt")
    output:
        families_to_species_strains = os.path.join(RUN_DIR,"database_resources/peptidotyping/families_to_species_strains.txt")
    log: os.path.join(RUN_DIR,"logs/peptidotyping/map_families_to_species_strains.log")
    container:
        config["containers"]["taxonkit"]
    shell:
        """
        export TAXONKIT_DB="{config[taxonkit_db_dir]}"

        # Fetch species and strains for each family taxid (skip header row).
        # If no families were detected, write an empty file and exit cleanly.
        ids=$(tail -n +2 {input.ncbi_taxonomy_ids} | awk '{{print $1}}')
        if [ -z "$ids" ]; then
            echo "No families detected; writing empty species/strain list." > {log}
            touch {output.families_to_species_strains}
        else
            echo "$ids" \
                | taxonkit list \
                | taxonkit filter -E species -E strain -N \
                > {output.families_to_species_strains} 2> {log}
        fi
        """

################################################################################
# Generate the second pass FASTA database
################################################################################
# Filter species_strain_lca_filtered_peptides.tsv to only include taxa that are
# descendants of the families detected in the first pass. This reduces the search
# space from all species/strains in UniProt to only those plausibly present.
rule generate_second_pass_db:
    input:
        species_tsv                 = os.path.join(config["peptidotyping_resource_dir"],"species_strain_lca_filtered_peptides.tsv"),
        families_to_species_strains = os.path.join(RUN_DIR,"database_resources/peptidotyping/families_to_species_strains.txt")
    output:
        second_pass_fasta = os.path.join(RUN_DIR,"database_resources/peptidotyping/second_pass_database.fasta")
    log: os.path.join(RUN_DIR,"logs/peptidotyping/generate_second_pass_db.log")
    container: config["containers"]["taxonkit"]
    shell:
        r"""
        set -euo pipefail
        log_ts() {{ echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a {log}; }}

        log_ts "Building taxid allowlist from families_to_species_strains.txt"

        ALLOWLIST=$(mktemp)
        awk 'NF && /^[[:space:]]*[0-9]/' {input.families_to_species_strains} \
            | tr -d '[:blank:]' | sort -u > "$ALLOWLIST"
        log_ts "Allowlist: $(wc -l < "$ALLOWLIST") species/strain taxids"

        log_ts "Filtering species_strain TSV to allowlist and writing FASTA"

        # TSV cols (1-indexed): id(1) sequence(2) lca(3) lca_il(4) fa(5) fa_il(6)
        #                       name(7) rank(8) parent_id(9) fasta_header(10)
        awk -F'\t' '
            ARGIND==1 {{ allowed[$1]=1; next }}
            FNR==1    {{ next }}
            $4 in allowed {{ print ">" $10; print $2 }}
        ' "$ALLOWLIST" {input.species_tsv} > {output.second_pass_fasta}

        rm -f "$ALLOWLIST"
        log_ts "Second-pass FASTA: $(grep -c "^>" {output.second_pass_fasta} || true) entries"
        """

################################################################################
# Performing the Second Pass Search, Resolve Species/Strains
################################################################################
# Search the filtered species/strain database with DIA-NN InfiniDIA to resolve
# which specific species and strains are present within the detected families.
rule perform_second_pass_search:
    input:
        raw_files_dir = os.path.join(EXPERIMENT_DIR,"input/ms_files"),
        fasta         = os.path.join(RUN_DIR,"database_resources/peptidotyping/second_pass_database.fasta"),
        config_file   = "config/peptidotyping_infinidia.cfg"
    output:
        second_pass_diann_parquet             = os.path.join(RUN_DIR,"database_resources/peptidotyping/second_pass_diann.parquet"),
        second_pass_diann_protein_description = os.path.join(RUN_DIR,"database_resources/peptidotyping/second_pass_diann.protein_description.tsv")
    log: os.path.join(RUN_DIR,"logs/peptidotyping/perform_second_pass_search.log")
    container: config["containers"]["diann"]
    threads: workflow.cores
    shell:
        """
        diann --cfg {input.config_file} \
            --fasta {input.fasta} \
            --out  {RUN_DIR}/database_resources/peptidotyping/second_pass_diann \
            --dir  {input.raw_files_dir} \
            --threads {threads} --verbose 1 >> {log} 2>&1
        """

################################################################################
# Infer the presence of species and strains based on the second pass search
# results using target-decoy FDR at the species/strain level.
################################################################################
rule infer_second_pass_presence:
    input:
        second_pass_diann = os.path.join(RUN_DIR,"database_resources/peptidotyping/second_pass_diann.parquet")
    output:
        detected_species_strains = os.path.join(RUN_DIR,"database_resources/peptidotyping/detected_species_strain_taxa_ids.txt"),
        fdr_results              = os.path.join(RUN_DIR,"database_resources/peptidotyping/second_pass_fdr_results.tsv")
    log: os.path.join(RUN_DIR,"logs/peptidotyping/infer_species_strain_presence.log")
    container: config["containers"]["conduitr"]
    script: "scripts/infer_species_strain_presence.R"

################################################################################
# Generate ncbi_taxa_ids.txt for handoff to the ncbi_taxonomy_id workflow
################################################################################
# The second-pass inference produces a clean TSV of detected species/strain
# taxids. Reformat it into the ncbi_taxa_ids.txt expected by the ncbi_taxonomy
# module (single column: ncbi_taxonomy_id).
rule generate_peptidotyping_ncbi_taxa_ids:
    input:
        detected_species_strains = os.path.join(RUN_DIR,"database_resources/peptidotyping/detected_species_strain_taxa_ids.txt")
    output:
        ncbi_taxa_ids = os.path.join(RUN_DIR,"ncbi_taxa_ids.txt")
    log: os.path.join(RUN_DIR,"logs/peptidotyping/generate_peptidotyping_ncbi_taxa_ids.log")
    shell:
        r"""
        set -euo pipefail
        echo "ncbi_taxonomy_id" > {output.ncbi_taxa_ids}
        tail -n +2 {input.detected_species_strains} \
            | cut -f1 | sort -u \
            >> {output.ncbi_taxa_ids} 2> {log}
        """

# ncbi_taxa_ids.txt is now consumed by the ncbi_taxonomy_id workflow.

