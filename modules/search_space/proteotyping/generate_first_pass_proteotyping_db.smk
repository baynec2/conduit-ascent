import glob
import os
EXPERIMENT_DIR = os.path.join("experiments",config["experiment"])
RAW_FILEPATHS = glob.glob(os.path.join(EXPERIMENT_DIR, "input/raw_files/*.raw"))
################################################################################
# Generating the Sequence Index
################################################################################
# This is needed to generate the file containing all peptides in TREMBL and SWISSPROT
# And their LCAs. See https://github.com/unipept/unipept-database/issues/75
rule build_sequence_index:
    output:
        "resources/databases/sequences.tsv.lz4",
        "resources/databases/relnotes.txt"
    container: "docker://baynec2/umgap:alpha"
    log: "resources/databases/logs/build_sequence_index.log"
    shell:
        """
        # Download the relnotes.txt file corresponding to the version of uniprotkb used
        curl -o resources/databases/relnotes.txt https://ftp.uniprot.org/pub/databases/uniprot/relnotes.txt 

        modules/search_space/proteotyping/scripts/unipept-database/scripts/generate_umgap_tables.sh tryptic \
            --output-dir resources/databases \
            --database-sources 'swissprot,trembl' \
            --temp-dir resources/databases/temp/ \
            --min-peptide-length 5 \
            --max-peptide-length 50 \
            >> {log} 2>&1
        """
################################################################################
# Determining Version of Uniprotkb that is being used in experiment
################################################################################
# The user should be able to know what version of uniprotkb was used to generate
# the sequence index. Furthermore, the user should be identified if the version
# of the sequence index is up to date with the current version of uniprotkb or not. 
rule check_sequence_index_version:
    input:
        "resources/databases/relnotes.txt"
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
# Generating the Proteotyping Fasta DB
################################################################################
# From the file containing the peptides and their LCAS, we can generate a fasta file
# only containing proteotypic peptides for a given taxonomy. We can subsequently use 
# this as a first pass database to identify taxa that are likely present in the experiment.
rule generate_first_pass_proteotyping_db:
    input:
        sequences = "resources/databases/sequences.tsv.lz4",
        taxons = "resources/databases/taxons.tsv.lz4"
    output:
        lca_filtered_taxa = "resources/databases/lca_filtered_taxa.tsv",
        first_pass_fasta = "resources/databases/first_pass_database.fasta"
    container: "docker://baynec2/conduitr:alpha"
    log: "resources/databases/logs/generate_first_pass_proteotyping_db.log"
    script: 
        "scripts/generate_first_pass_proteotyping_db.R"
################################################################################
# Generating the First Pass Spectral Library
################################################################################
# We can generate a spectral library from the proteotyping fasta db to use with
# Diann. This can then be used to perform the first pass search identifying 
# species that are likely present in the experiment. 
rule generate_first_pass_spectral_library:
    input: 
        fasta = "resources/databases/first_pass_database.fasta",
        config_file = "config/proteotyping_firstpass_diann_spectral_library.cfg"
    output: "resources/databases/first_pass_database.predicted.speclib"
    container: "docker://baynec2/diann2.1.0:alpha"
    log: "resources/databases/logs/generate_first_pass_spectral_library.log"
    threads: workflow.cores
    shell:
        """
        diann --cfg {input.config_file} \
        --fasta {input.fasta} \
        --threads {threads} \
        --out-lib resources/databases/first_pass_database >> {log} 2>&1
        """
################################################################################
# Performing the First Pass Search
################################################################################
# Searching our first pass spectral library with Diann
rule perform_first_pass_search:
    input:
        raw_files_dir = os.path.join(EXPERIMENT_DIR,"input/raw_files"),
        spectral_library = "resources/databases/first_pass_database.predicted.speclib",
        fasta = "resources/databases/first_pass_database.fasta",
        config_file = "config/proteotyping_firstpass_diann.cfg"
    output:
        first_pass_diann_parquet = os.path.join(EXPERIMENT_DIR,"input/database_resources/proteotyping/first_pass_diann.parquet"),
        first_pass_diann_protein_description =  os.path.join(EXPERIMENT_DIR,"input/database_resources/proteotyping/first_pass_diann.protein_description.tsv")
    log: os.path.join(EXPERIMENT_DIR,"logs/proteotyping/perfrom_first_pass_search.log")
    container:
        "docker://baynec2/diann2.1.0:alpha"
    threads: workflow.cores 
    shell:
        """
        diann --cfg {input.config_file} \
        --fasta {input.fasta} \
        --out  experiments/{config[experiment]}/input/database_resources/proteotyping/first_pass_diann \
        --dir {input.raw_files_dir} \
        --lib {input.spectral_library} \
        --threads {threads} --verbose 1 >> {log} 2>&1
        """

################################################################################
# Scoring the likelihood of species being present based on the first pass search results
################################################################################
# We aim to score the likelihood that a species is truly present based on the number of detected
# species-specific proteotypic peptides. To do this, we calculate the fraction of that species’ 
# known proteotypic peptides that were observed in the experiment, and compare this to the number 
# we would expect to appear by chance under a 1% FDR.
# The current implentation is very basic, and the statistical framework needs to be improved. 
rule infer_species_presence:
    input:
        first_pass_diann = os.path.join(EXPERIMENT_DIR,"input/database_resources/proteotyping/first_pass_diann.parquet"),
        taxon_specific_peptide_db =  "resources/databases/lca_filtered_taxa.tsv"
    output:
        first_pass_search_taxa_metrics = os.path.join(EXPERIMENT_DIR,"input/database_resources/proteotyping/first_pass_search_taxa_metrics.tsv"),
        called_taxa_per_run = os.path.join(EXPERIMENT_DIR,"input/database_resources/proteotyping/called_taxonomy_per_run.tsv"),
        ncbi_taxonomy_id = os.path.join(EXPERIMENT_DIR,"input/ncbi_taxa_ids.txt")
    log: os.path.join(EXPERIMENT_DIR,"logs/proteotyping/infer_species_presence.log")
    script: "scripts/infer_species_presence.R"

# After this, the ncbi_taxa_ids.txt file will get plugged into the ncbi_taxonomy_id workflow and the second pass search will start.


   