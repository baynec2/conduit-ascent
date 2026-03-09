import glob
import os
EXPERIMENT_DIR = os.path.join("experiments",config["experiment"])
RAW_FILEPATHS = glob.glob(os.path.join(EXPERIMENT_DIR, "input/raw_files/*.raw"))

################################################################################
# Generating the Sequence Index
################################################################################
# This is needed to generate the file containing all peptides in TREMBL and SWISSPROT
# And their LCAS. See https://github.com/unipept/unipept-database/issues/75
rule build_sequence_index:
    output:
        sequences = os.path.join(config["peptidotyping_resource_dir"],"sequences.tsv.lz4"),
        taxons    = os.path.join(config["peptidotyping_resource_dir"],"taxons.tsv.lz4")
    params:
        outdir = config["peptidotyping_resource_dir"],
        temp_outdir = os.path.join(config["peptidotyping_resource_dir"],"temp")
    log:
        os.path.join(config["peptidotyping_resource_dir"],"logs/build_sequence_index.log")
    container:
        "docker://baynec2/umgap:alpha"
    shell:
        r"""
        set -euo pipefail

        # Ensure directories exist
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})
        mkdir -p {params.temp_outdir}

        # Download UniProt release notes
        curl -L \
          -o {params.outdir}/relnotes.txt \
          https://ftp.uniprot.org/pub/databases/uniprot/relnotes.txt

        # Setting temp dir
        export TMPDIR={params.temp_outdir}

        # Build UMGAP peptidotyping tables
        modules/search_space/peptidotyping/scripts/unipept-database/scripts/generate_umgap_tables.sh tryptic \
          --output-dir {params.outdir} \
          --database-sources swissprot,trembl \
          --temp-dir {params.temp_outdir} \
          --min-peptide-length 5 \
          --max-peptide-length 50 \
          >> {log} 2>&1
        """
################################################################################
# Determining Version of Uniprotkb that is being used in experiment
################################################################################
rule check_sequence_index_version:
    input:
        relnotes = os.path.join(config["peptidotyping_resource_dir"],"relnotes.txt")
    output:
        version_file = os.path.join(config["peptidotyping_resource_dir"],".version")
    params:
        outdir = config["peptidotyping_resource_dir"]
    log:
        os.path.join(EXPERIMENT_DIR,"logs/search_space/check_sequence_index_version.log")
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
    container: "docker://baynec2/conduitr:alpha"
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
    container: "docker://baynec2/conduitr:alpha"
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
# Determining which families lack family-specific peptides
################################################################################
# Some families do not have peptides that are specific at the family level
# when using LCA-based assignment. An empirically observed example is
# Akkermansiaceae, which is monophyletic and contains a single genus.
# In this case, all taxon-informative peptides resolve to the genus level
# by definition. As a result, a strategy that relies solely on family-level
# LCA peptides would fail to detect descendant species within such clades.
#
# To address this, we first identify families that lack family-specific
# peptides. For these cases, we generate a mapping to peptides whose LCA
# is at the genus level, allowing the workflow to fall back to genus-level
# assignments when appropriate.

# First we need to determine all the families that we can detect (with a strain/species corresponding)
rule determine_all_possible_families:
    input:
        species_strain_lca_filtered_peptides = os.path.join(config["peptidotyping_resource_dir"],"species_strain_lca_filtered_peptides.tsv")
    output:
        possible_family_taxons = os.path.join(config["peptidotyping_resource_dir"],"possible_family_taxons.txt")
    threads: 8
    container: "docker://baynec2/conduitr:alpha"
    shell: """
        tail -n +2 {input.species_strain_lca_filtered_peptides} | cut -f4 | awk 'NF && $1!=""' | sort -u \
        | taxonkit lineage -j {threads} \
        | taxonkit reformat -t -r -f '{{taxid}}\t{{ranks}}\t{{lineage_taxids}}' \
        | awk -F'\t' 'NR>1 {{ n=split($2,rk,";"); split($3,ln,";"); for(i=1;i<=n;i++) if(rk[i]=="family" && ln[i]!="") print ln[i] }}' \
        | sort -u > {output.possible_family_taxons}
        """
# Next we need to determine what families were actually detected.
rule determine_detected_families:
    input:
        family_lca_filtered_peptides = os.path.join(config["peptidotyping_resource_dir"],"family_lca_filtered_peptides.tsv"),
    output:
        detected_families = os.path.join(config["peptidotyping_resource_dir"],"detected_families.txt")
    container: "docker://baynec2/conduitr:alpha"
    shell: """
        tail -n +2 {input.family_lca_filtered_peptides} | cut -f4 | awk 'NF && $1!=""' | sort -u > {output.detected_families}
        """
# Now we need to determine what families were missing (in possible but not in detected).
rule determine_missing_families:
    input:
        possible_family_taxons = os.path.join(config["peptidotyping_resource_dir"],"possible_family_taxons.txt"),
        detected_families = os.path.join(config["peptidotyping_resource_dir"],"detected_families.txt"),
    output:
        missing_families = os.path.join(config["peptidotyping_resource_dir"],"missing_families.txt")
    container: "docker://baynec2/conduitr:alpha"
    shell: """
        comm -23 <(sort {input.possible_family_taxons}) <(sort {input.detected_families}) > {output.missing_families}
        """
# For each missing family, find genus-level peptides that belong to that family (fallback for detection).
rule find_fallback_genus_for_missing_families:
    input:
        missing_families = os.path.join(config["peptidotyping_resource_dir"],"missing_families.txt"),
        species_strain_lca_filtered_peptides = os.path.join(config["peptidotyping_resource_dir"],"species_strain_lca_filtered_peptides.tsv"),
    output:
        missing_families_genus_fallback = os.path.join(config["peptidotyping_resource_dir"],"missing_families_genus_fallback.tsv")
    threads: 8
    container: "docker://baynec2/conduitr:alpha"
    shell: """
        set -e
        # Build lca_il -> rank,family_taxid,genus_taxid for all unique lca_il in species/strain peptides
        tail -n +2 {input.species_strain_lca_filtered_peptides} | cut -f4 | awk 'NF && $1!=""' | sort -u \
        | taxonkit lineage -j {threads} \
        | taxonkit reformat -t -r -f '{{taxid}}\t{{rank}}\t{{ranks}}\t{{lineage_taxids}}' \
        | awk -F'\t' 'NR>1 {{
          n=split($3,rk,";"); split($4,ln,";");
          r=$2; fam=""; gen="";
          for(i=1;i<=n;i++) {{
            if(rk[i]=="family") fam=ln[i];
            if(rk[i]=="genus")  gen=ln[i];
          }}
          if(fam!="" && gen!="") print $1"\t"r"\t"fam"\t"gen
        }}' > lca_il_to_fam_genus.tsv
        # Join: output family_taxid, genus_taxid, peptide_id for peptides with LCA rank genus and family in missing_families
        echo -e "family_taxid\tgenus_taxid\tpeptide_id" > {output.missing_families_genus_fallback}
        awk -F'\t' -v OFS='\t' '
          NR==FNR {{ miss[$1]=1; next }}
          FILENAME!=prev {{ prev=FILENAME; f++ }}
          f==1 {{ rank[$1]=$2; fam[$1]=$3; gen[$1]=$4; next }}
          f==2 && FNR==1 {{ next }}
          f==2 {{ lca=$4; id=$1; if(lca in rank && rank[lca]=="genus" && fam[lca] in miss) print fam[lca], gen[lca], id }}
        ' {input.missing_families} lca_il_to_fam_genus.tsv {input.species_strain_lca_filtered_peptides} >> {output.missing_families_genus_fallback}
        """    

################################################################################
# Generating the First Pass Spectral Library
################################################################################
# We can generate a spectral library from the peptidotyping fasta db to use with
# Diann. This can then be used to perform the first pass search identifying 
# species that are likely present in the experiment. 
rule generate_first_peptidotyping_spectral_library:
    input: 
        fasta = os.path.join(config["peptidotyping_resource_dir"],"phylum_peptidotyping_db.fasta"),
        config_file = "config/peptidotyping_firstpass_diann_spectral_library.cfg"
    output: 
        os.path.join(config["peptidotyping_resource_dir"],"phylum_peptidotyping.predicted.speclib")
    # DIANN adds the .predicted.speclib extennsion 
    container: "docker://baynec2/diann2.1.0:alpha"
    log: os.path.join(config["peptidotyping_resource_dir"],"logs/generate_phylumn_peptidotyping_spectral_library.log")
    threads: workflow.cores
    shell:
        """
        diann --cfg {input.config_file} \
        --fasta {input.fasta} \
        --threads {threads} \
        --out-lib {config[peptidotyping_resource_dir]}phylum_peptidotyping >> {log} 2>&1
        """
################################################################################
# Generating a first pass database of species specific (or more granular)- likley
# highly abundant proteins based on GO term functional annotation
################################################################################
# Here we will generate a first pass database consisting of peptides with a LCA
# at the species level or higher that is constrained to only include GO Terms that 
# are likely to be highly abundant, and thus detectable. 
# This should work because if highly abundant protiens are not detected, much less
# abundant proteins are extremely unlikely to be detected as well.
# We will use the following GO Terms to restrict the database 
## GO TERMS
#GO:0005840    ribosome 22,199,805 annotations
#GO:0006412    translation 23,081,373 annotations
#GO:0003746    translation elongation factor activity 437,747 annotations
#GO:0005856    cytoskeleton 4,967,930 annotations
#GO:0008152    metabolic process 139,512,380 annotations
#GO:0016020    membrane 62,787,788 annotations
#GO:0003677    DNA binding 21,017,894 annotations
#GO:0003723    RNA binding 18,400,347 annotations
# This is similar conceptually to the HAPiID approach, but is applied to a more general set of peptides.
# in our implementation, see https://pmc.ncbi.nlm.nih.gov/articles/PMC8017886/ for HAPiID paper.

# Or alternatively:
#GO:0005840  # ribosome
#GO:0006412  # translation
#GO:0003746  # translation elongation factor activity
#GO:0006457  # protein folding (chaperones)
#GO:0051082  # unfolded protein binding (chaperones)
#GO:0016887  # ATPase activity
#GO:0006260  # DNA replication
#GO:0003677  # DNA binding
#GO:0003723  # RNA binding
#GO:0016020  # membrane
#GO:0005198  # structural molecule activity
#GO:0005856  # cytoskeleton 4,967,930 annotations

rule generate_highly_abundant_peptidotyping_database:
    input:
        sequences = os.path.join(config["peptidotyping_resource_dir"],"sequences.tsv.lz4"),
        taxons = os.path.join(config["peptidotyping_resource_dir"],"taxons.tsv.lz4")
    output:
        lca_filtered_taxa = os.path.join(config["peptidotyping_resource_dir"],"translation_lca_filtered_peptides.tsv"),
        first_pass_fasta = os.path.join(config["peptidotyping_resource_dir"],"translation_peptidotyping_db.fasta")
    params:
        taxon_ranks_str = "species,strain",  # Comma-separated string of ranks
        go_terms = "GO:0005840,GO:0006412,GO:0003746"  # Comma-separated GO terms
    container: "docker://baynec2/conduitr:alpha"
    log: os.path.join(config["peptidotyping_resource_dir"],"logs/generate_translation_peptidotyping_database.log")
    shell:
        r"""
        set -euo pipefail
        
        SEQUENCES_FILE="{input.sequences}"
        TAXONS_FILE="{input.taxons}"
        LCA_FILTERED_TAXA="{output.lca_filtered_taxa}"
        OUTPUT_FASTA="{output.first_pass_fasta}"
        TAXON_RANKS="{params.taxon_ranks_str}"
        GO_TERMS="{params.go_terms}"
        LOG_FILE="{log}"
        
        # Function to log with timestamp
        log_with_timestamp() {{
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
        }}
        
        log_with_timestamp "Starting generate_highly_abundant_peptidotyping_database"
        log_with_timestamp "Input sequences: $SEQUENCES_FILE"
        log_with_timestamp "Input taxons: $TAXONS_FILE"
        log_with_timestamp "Output TSV: $LCA_FILTERED_TAXA"
        log_with_timestamp "Output FASTA: $OUTPUT_FASTA"
        log_with_timestamp "Taxon ranks: $TAXON_RANKS"
        log_with_timestamp "GO terms: $GO_TERMS"
        
        # Convert comma-separated ranks to awk pattern
        RANK_PATTERN=$(echo "$TAXON_RANKS" | tr ',' '|')
        
        # Convert comma-separated GO terms to awk pattern (for matching in JSON)
        # We'll search for patterns like "GO:0005840" in the JSON string
        # Format: "GO:0005840|GO:0006412|..." - matches any of the GO terms
        GO_PATTERN=$(echo "$GO_TERMS" | tr ',' '|')
        
        # Create GO terms array file for tracking individual GO terms
        echo "$GO_TERMS" | tr ',' '\n' > /tmp/go_terms_list_$$.txt
        
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
        
        log_with_timestamp "Step 2: Streaming sequences, filtering by taxonomy AND GO terms, writing output (single pass)"
        
        # Create GO terms list file for tracking
        echo "$GO_TERMS" | tr ',' '\n' > /tmp/go_terms_list_$$.txt
        
        # Single pass - stream sequences, filter by taxonomy AND GO terms, write TSV + FASTA
        awk -F'\t' \
            -v rank_pattern="$RANK_PATTERN" \
            -v go_pattern="$GO_PATTERN" \
            -v go_terms_file="/tmp/go_terms_list_$$.txt" \
            -v tsv_file="$LCA_FILTERED_TAXA" \
            -v fasta_file="$OUTPUT_FASTA" \
            -v log_file="$LOG_FILE" \
            '
            BEGIN {{
                # Load individual GO terms for tracking statistics
                go_term_count = 0
                while ((getline go_term < go_terms_file) > 0) {{
                    go_term_count++
                    go_terms_list[go_term_count] = go_term
                    go_term_peptide_count[go_term] = 0
                }}
                close(go_terms_file)
            }}
            
            # Process first file (taxons): build lookup
            FNR == NR {{
                if ($3 ~ "^(" rank_pattern ")$") {{
                    taxon_name[$1] = $2
                    taxon_rank[$1] = $3
                    taxon_parent[$1] = $4
                }}
                next
            }}
            
            # Process second file (sequences): filter by taxonomy AND GO terms, then write
            {{
                total++
                
                seq_id = $1
                sequence = $2
                lca = $3
                lca_il = $4
                fa = $5
                fa_il = $6  # This column contains JSON with GO terms
                
                # Filter 1: Check if lca_il matches taxonomy
                taxonomy_match = (lca_il in taxon_name)
                
                # Filter 2: Check if fa_il (JSON) contains any of the specified GO terms
                # Also track which specific GO terms matched for statistics
                # The JSON format is: {{"num":{{...}},"data":{{"GO:0005840":2,"GO:0006412":1,...}}}}
                go_match = 0
                if (fa_il != "" && fa_il != "\\N" && taxonomy_match) {{
                    # Check each GO term individually to track per-GO-term statistics
                    for (i = 1; i <= go_term_count; i++) {{
                        go_term = go_terms_list[i]
                        go_regex = "\"" go_term "\""
                        if (match(fa_il, go_regex)) {{
                            go_match = 1
                            go_term_peptide_count[go_term]++
                            # Track unique taxa per GO term (using composite key)
                            go_term_taxa_key = go_term SUBSEP lca_il
                            go_term_taxa[go_term_taxa_key] = 1
                        }}
                    }}
                }}
                
                # Keep only if BOTH filters pass
                if (taxonomy_match && go_match) {{
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
                printf "[%s] Filtered to %d peptides with LCA at specified ranks AND matching GO terms (from %d total, %.2f%%)\n", strftime("%Y-%m-%d %H:%M:%S"), kept, total, (kept/total)*100 >> log_file
                
                # Output GO term statistics
                printf "\n[%s] GO Term Statistics:\n", strftime("%Y-%m-%d %H:%M:%S") >> log_file
                printf "GO_Term\tPeptides\tUnique_Taxa\n" >> log_file
                for (i = 1; i <= go_term_count; i++) {{
                    go_term = go_terms_list[i]
                    # Count unique taxa for this GO term
                    taxa_count = 0
                    for (key in go_term_taxa) {{
                        split(key, arr, SUBSEP)
                        if (arr[1] == go_term) {{
                            taxa_count++
                        }}
                    }}
                    printf "%s\t%d\t%d\n", go_term, go_term_peptide_count[go_term], taxa_count >> log_file
                }}
                close(log_file)
            }}
            ' <(lz4 -d -c "$TAXONS_FILE") <(lz4 -d -c "$SEQUENCES_FILE")
        
        # Cleanup
        rm -f /tmp/go_terms_list_$$.txt
        
        log_with_timestamp "Completed generate_highly_abundant_peptidotyping_database"
        log_with_timestamp "TSV file: $LCA_FILTERED_TAXA"
        log_with_timestamp "FASTA file: $OUTPUT_FASTA"
        """
################################################################################
# Performing the First Pass Search
################################################################################
# Searching our first pass spectral library with Diann
rule perform_first_pass_search:
    input:
        raw_files_dir = os.path.join(EXPERIMENT_DIR,"input/raw_files"),
        spectral_library = os.path.join(config["peptidotyping_resource_dir"],"first_pass_database.predicted.speclib"),
        fasta = os.path.join(config["peptidotyping_resource_dir"],"first_pass_database.fasta"),
        config_file = "config/proteotyping_firstpass_diann.cfg"
    output:
        first_pass_diann_parquet = os.path.join(EXPERIMENT_DIR,"input/database_resources/peptidotyping/first_pass_diann.parquet"),
        first_pass_diann_protein_description =  os.path.join(EXPERIMENT_DIR,"input/database_resources/peptidotyping/first_pass_diann.protein_description.tsv")
    log: os.path.join(EXPERIMENT_DIR,"logs/peptidotyping/perfrom_first_pass_search.log")
    container:
        "docker://baynec2/diann2.1.0:alpha"
    threads: workflow.cores 
    shell:
        """
        diann --cfg {input.config_file} \
        --fasta {input.fasta} \
        --out  experiments/{config[experiment]}/input/database_resources/peptidotyping/first_pass_diann \
        --dir {input.raw_files_dir} \
        --lib {input.spectral_library} \
        --threads {threads} --verbose 1 >> {log} 2>&1
        """
################################################################################
# Determine what families are present based on the first pass results
################################################################################
# Here we will just use the number of peptides that were detected to infer the presence of families. 
# If a number of peptides > than the theshold are found, we will infer that the family is present and 
# only use the species/strain level peptides that belong to the family. 
rule infer_family_presence:
    input:
        first_pass_diann = os.path.join(EXPERIMENT_DIR,"input/database_resources/peptidotyping/first_pass_diann.parquet"),
    output:
        ncbi_taxonomy_id = os.path.join(EXPERIMENT_DIR,"input/detected_family_taxa_ids.txt")
    params:
        threshold = 2
    log: os.path.join(EXPERIMENT_DIR,"logs/peptidotyping/infer_family_presence.log")
    script: "scripts/infer_family_presence.R"

################################################################################
# Determine which species and strains belong to the detected families. 
################################################################################
# By mapping these values, we can tell what species there is family level evidence for.
# This will allow us to subsequently search a reduced strain/species specific peptide database. 
 rule map_families_to_species_strains:
    input:
        ncbi_taxonomy_ids = os.path.join(EXPERIMENT_DIR,"input/database_resources/peptidotyping/detected_family_taxa_ids.txt")
    output:
        families_to_species_strains = os.path.join(EXPERIMENT_DIR,"input/database_resources/peptidotyping/families_to_species_strains.txt")
    log: os.path.join(EXPERIMENT_DIR,"logs/peptidotyping/map_families_to_species_strains.log")
    container:
        "quay.io/biocontainers/taxonkit:0.20.0--h9ee0642_1"
    shell:
        """
        # Make sure Taxonkit knows where the local taxonomy DB (optional)
        export TAXONKIT_DB=${TAXONKIT_DB:-/root/.taxonkit}  # or mount your local DB if needed

        # Fetch species and strains for each family taxid
        cut -f1 {input.ncbi_taxonomy_ids} | xargs -I {{}} taxonkit list --id {{}} --rank species --rank strain > {output.families_to_species_strains} 2> {log}
        """


# Now we need to filter the peptidotyping database to only include species and strain level peptides from the families that were detected. 
rule filter_species_strain_peptidotyping_database:
    input:
        species_strain_fasta = os.path.join(config["peptidotyping_resource_dir"],"species_strain_peptidotyping_db.fasta")
    output:
        filtered_species_strain_fasta = os.path.join(config["peptidotyping_resource_dir"],"filtered_species_strain_peptidotyping_db.fasta")
    log: os.path.join(config["peptidotyping_resource_dir"],"logs/filter_species_strain_peptidotyping_database.log")
    shell:
        """
        grep -f {input.species_strain_fasta} {input.species_strain_fasta} > {output.filtered_species_strain_fasta}
        """

# After this, the ncbi_taxa_ids.txt file will get plugged into the ncbi_taxonomy_id workflow and the second pass search will start.

