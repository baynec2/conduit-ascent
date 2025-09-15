EXPERIMENT_DIR = os.path.join("experiments",config["experiment"])
################################################################################
# Getting Database Resources (Defining the Search Space)
################################################################################        
rule get_fasta:
    input:
        os.path.join(EXPERIMENT_DIR,"input/ncbi_taxa_ids.txt")
    output:
        # File containing proteome ids. Allows user to see what proteomes for each organism were downloaded.
        # Also contains NA values for NCBI taxon ids that are not in the Uniprot database.
       os.path.join(EXPERIMENT_DIR,"input/database_resources/proteome_ids.txt"),
       os.path.join(EXPERIMENT_DIR,"input/database_resources/database.fasta")
    log: os.path.join(EXPERIMENT_DIR,"logs/search_space/ncbi_taxonomy/get_fasta.log")
    container:
        "docker://baynec2/conduitr:alpha"
    script:
        "scripts/get_fasta.R"

# Generating Full Taxonomy Information for each Organism ID
rule get_taxonomy:
    input:
        os.path.join(EXPERIMENT_DIR,"input/ncbi_taxa_ids.txt"),
        # File containing results of get fasta. Includes info about what taxa ids had cooresponding proteome ids.
        os.path.join(EXPERIMENT_DIR,"input/database_resources/proteome_ids.txt")
    output:
        os.path.join(EXPERIMENT_DIR,"input/database_resources/taxonomy.txt")
    log: os.path.join(EXPERIMENT_DIR,"logs/search_space/ncbi_taxonomy/get_taxonomy.log")
    container: "docker://baynec2/conduitr:alpha"
    script:
      "scripts/get_taxonomy.R"
      
# Taking the data from fasta and taxonomy file, putting it in fasta.
rule get_protein_info_from_fasta:
    input:
      database_fasta= os.path.join(EXPERIMENT_DIR,"input/database_resources/database.fasta"),
      taxonomy_txt= os.path.join(EXPERIMENT_DIR,"input/database_resources/taxonomy.txt")
    output:
      os.path.join(EXPERIMENT_DIR,"input/database_resources/protein_info.txt")
    log:
      os.path.join(EXPERIMENT_DIR,"logs/search_space/ncbi_taxonomy/get_protein_info_from_fasta.log")
    container: "docker://baynec2/conduitr:alpha"
    script:
      "scripts/get_protein_info_from_fasta.R"
      
# Plotting a taxonomic tree containing the taxonomy used in experiment
rule plot_taxonomic_tree:
    input: os.path.join(EXPERIMENT_DIR,"input/database_resources/taxonomy.txt")
    output: os.path.join(EXPERIMENT_DIR,"input/database_resources/taxonomic_tree_of_database.pdf")
    log: os.path.join(EXPERIMENT_DIR,"logs/search_space/ncbi_taxonomy/plot_taxonomic_tree.log")
    container:"docker://baynec2/conduitr:alpha"
    script:
      "scripts/plot_taxonomic_tree.R"

# Creating a DataBase ReadMe          
rule make_database_resources_readme:
    input: os.path.join(EXPERIMENT_DIR,"input/database_resources/protein_info.txt")
    output:
      md = os.path.join(EXPERIMENT_DIR,"input/database_resources/README.md"),
      html = os.path.join(EXPERIMENT_DIR,"input/database_resources/README.html")
    log: os.path.join(EXPERIMENT_DIR,"logs/search_space/ncbi_taxonomy/make_database_resources_readme.log")
    container: "docker://baynec2/conduitr:alpha"
    script:
     "scripts/make_database_resources_readme.R"