EXPERIMENT_DIR = os.path.join("experiments",config["experiment"])
################################################################################
# Processing Matrices
################################################################################
rule process_diann_matrices:
  input:
    report_pg_matrix=os.path.join(EXPERIMENT_DIR,"output/diann_output/diann.pg_matrix.tsv"),
    report_pr_matrix =os.path.join(EXPERIMENT_DIR,"output/diann_output/diann.pr_matrix.tsv")
  output: 
    protein_group_matrix = os.path.join(EXPERIMENT_DIR,"output/output_files/protein_group_matrix.tsv"),
    precursor_matrix = os.path.join(EXPERIMENT_DIR,"output/output_files/precursor_matrix.tsv"),
    peptide_matrix = os.path.join(EXPERIMENT_DIR,"output/output_files/peptide_matrix.tsv")
  log: os.path.join(EXPERIMENT_DIR,"logs/matrices/process_diann_matrices.log")
  container: "apptainer/conduitR.sif"
  script: "scripts/process_diann_matrices.R"
    
rule process_taxonomic_matrices:
  input:
    report_pr_matrix=os.path.join(EXPERIMENT_DIR,"output/output_files/precursor_matrix.tsv"),
    report_pg_matrix=os.path.join(EXPERIMENT_DIR,"output/output_files/protein_group_matrix.tsv"),
    protein_info=os.path.join(EXPERIMENT_DIR,"input/database_resources/protein_info.txt")
  output:
    domain_matrix=os.path.join(EXPERIMENT_DIR,"output/output_files/domain_matrix.tsv"),
    kingdom_matrix=os.path.join(EXPERIMENT_DIR,"output/output_files/kingdom_matrix.tsv"),
    phylum_matrix=os.path.join(EXPERIMENT_DIR,"output/output_files/phylum_matrix.tsv"),
    class_matrix=os.path.join(EXPERIMENT_DIR,"output/output_files/class_matrix.tsv"),
    order_matrix=os.path.join(EXPERIMENT_DIR,"output/output_files/order_matrix.tsv"),
    family_matrix=os.path.join(EXPERIMENT_DIR,"output/output_files/family_matrix.tsv"),
    genus_matrix=os.path.join(EXPERIMENT_DIR,"output/output_files/genus_matrix.tsv"),
    species_matrix=os.path.join(EXPERIMENT_DIR,"output/output_files/species_matrix.tsv")
  log: os.path.join(EXPERIMENT_DIR,"logs/matrices/process_taxonomic_matrices.log")
  container: "apptainer/conduitR.sif"
  script: "scripts/process_taxonomic_matrices.R"
    
rule process_go_matrix:
  input:
    report_pg_matrix=os.path.join(EXPERIMENT_DIR,"output/output_files/protein_group_matrix.tsv"),
    go_annotations = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/go_annotations.txt")
  output: 
    go_annotations_matrix = os.path.join(EXPERIMENT_DIR,"output/output_files/go_matrix.tsv"),
    go_annotations_taxa_matrix = os.path.join(EXPERIMENT_DIR,"output/output_files/go_taxa_matrix.tsv")
  log: os.path.join(EXPERIMENT_DIR,"logs/matrices/process_go_matrix.log")
  container: "apptainer/conduitR.sif"
  script: "scripts/process_go_matrix.R"
  
rule process_subcellular_locations_matrix:
  input:
    report_pg_matrix=os.path.join(EXPERIMENT_DIR,"output/output_files/protein_group_matrix.tsv"),
    subcellular_locations = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/subcellular_locations.txt")
  output: 
    subcellular_locations_matrix = os.path.join(EXPERIMENT_DIR,"output/output_files/subcellular_locations_matrix.tsv")
  log: os.path.join(EXPERIMENT_DIR,"logs/matrices/process_subcellular_locations_matrix.log")
  container: "apptainer/conduitR.sif"
  script:"scripts/process_subcellular_locations_matrix.R"
  
rule process_kegg_matrices:
  input:
      pg_matrix=os.path.join(EXPERIMENT_DIR,"output/output_files/protein_group_matrix.tsv"),
      kegg_annotations = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/kegg_annotations.txt")
  output: 
    kegg_pathway_matrix = os.path.join(EXPERIMENT_DIR,"output/output_files/kegg_pathway_matrix.tsv"),
    kegg_ko_matrix = os.path.join(EXPERIMENT_DIR,"output/output_files/kegg_ko_matrix.tsv")
  log: os.path.join(EXPERIMENT_DIR,"logs/matrices/process_kegg_matrices.log")
  container: "apptainer/conduitR.sif"
  script:"scripts/process_kegg_matrices.R"  
 

# Moving Database Resources to Output Directory.
rule move_database_resources:
    input:
        expand(os.path.join(EXPERIMENT_DIR,"input/database_resources/{filename}"),
               filename=[
                   "database.fasta",
                   "proteome_ids.txt",
                   "taxonomy.txt",
                   "protein_info.txt",
                   "taxonomic_tree_of_database.pdf",
                   "database.predicted.speclib",
                   "README.md",
                   "README.html"
               ]),
        expand(os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/{filename}"), 
               filename=[
                   "detected_protein_info.txt",
                   "detected_protein.fasta",
                   "uniprot_annotated_protein_info.txt",
                   "go_annotations.txt",
                   "subcellular_locations.txt",
                   "kegg_annotations.txt"
               ]),
    output:
        expand(os.path.join(EXPERIMENT_DIR,"output/database_resources/{filename}"), 
               filename=[
                   "database.fasta",
                   "proteome_ids.txt",
                   "taxonomy.txt",
                   "protein_info.txt",
                   "taxonomic_tree_of_database.pdf",
                   "database.predicted.speclib",
                   "README.md",
                   "README.html"
               ]),
        expand(os.path.join(EXPERIMENT_DIR,"output/database_resources/detected_protein_resources/{filename}"), 
               filename=[
                   "detected_protein_info.txt",
                   "detected_protein.fasta",
                   "uniprot_annotated_protein_info.txt",
                   "go_annotations.txt",
                   "subcellular_locations.txt",
                   "kegg_annotations.txt"
               ]),
    log: os.path.join(EXPERIMENT_DIR,"logs/matrices/move_database_resources.log")
    shell:
        """
        mkdir -p {EXPERIMENT_DIR}/output/database_resources/detected_protein_resources
        cp -u -r {EXPERIMENT_DIR}/input/database_resources/* {EXPERIMENT_DIR}/output/database_resources/
        cp -u -r {EXPERIMENT_DIR}/input/database_resources/detected_protein_resources/* {EXPERIMENT_DIR}/output/database_resources/detected_protein_resources/
        """