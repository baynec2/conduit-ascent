import os
import glob
# Experiment specific directories
EXPERIMENT_DIR = config["experiment_dir"]
RUN_DIR = config["run_dir"]
# Taking the data from fasta and taxonomy file, putting it in fasta.
rule get_protein_info_from_fasta:
    input:
      database_fasta= os.path.join(RUN_DIR,"database_resources/database.fasta"),
      taxonomy_txt= os.path.join(RUN_DIR,"database_resources/taxonomy.txt")
    output:
      os.path.join(RUN_DIR,"database_resources/protein_info.txt")
    log:
      os.path.join(RUN_DIR,"logs/search_space/database_processing/get_protein_info_from_fasta.log")
    container: "docker://baynec2/conduitr:alpha"
    script:
      "scripts/get_protein_info_from_fasta.R"

# Plotting a taxonomic tree containing the taxonomy used in experiment
rule plot_taxonomic_tree:
    input: os.path.join(RUN_DIR,"database_resources/taxonomy.txt")
    output: os.path.join(RUN_DIR,"database_resources/taxonomic_tree_of_database.pdf")
    log: os.path.join(RUN_DIR,"logs/search_space/database_processing/plot_taxonomic_tree.log")
    container:"docker://baynec2/conduitr:alpha"
    script:
      "scripts/plot_taxonomic_tree.R"

# Creating a DataBase ReadMe
rule make_database_resources_readme:
    input: os.path.join(RUN_DIR,"database_resources/protein_info.txt")
    output:
      md = os.path.join(RUN_DIR,"database_resources/README.md"),
      html = os.path.join(RUN_DIR,"database_resources/README.html")
    log: os.path.join(RUN_DIR,"logs/search_space/database_processing/make_database_resources_readme.log")
    container: "docker://baynec2/conduitr:alpha"
    script:
     "scripts/make_database_resources_readme.R"
