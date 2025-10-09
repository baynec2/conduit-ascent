EXPERIMENT_DIR = os.path.join("experiments",config["experiment"])
################################################################################
# Processing to QFeatures
################################################################################
rule generate_qfeatures_from_diann_parquet:
  input:
    diann_parquet=os.path.join(EXPERIMENT_DIR,"output/diann_output/diann.parquet"),
    sample_annotation=os.path.join(EXPERIMENT_DIR,"input/sample_annotation.txt")
  output: 
    qf = os.path.join(EXPERIMENT_DIR,"output/output_files/qf.rds")
  log: os.path.join(EXPERIMENT_DIR,"logs/build_conduit/generate_qfeatures_from_diann_parquet.log")
  container: "docker://baynec2/conduitr:alpha"
  script: "scripts/generate_qfeatures_from_diann_parquet.R"
    
rule add_annotations_to_qfeatures:
  input:
    qf = os.path.join(EXPERIMENT_DIR,"output/output_files/qf.rds"),
    uniprot_annotated_protein_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/uniprot_annotated_protein_info.txt")
  output:
    annotated_qf=os.path.join(EXPERIMENT_DIR,"output/output_files/annotated_qf.rds")
  log: os.path.join(EXPERIMENT_DIR,"logs/build_conduit/add_annotations_to_qfeatures.log")
  container: "docker://baynec2/conduitr:alpha"
  script: "scripts/add_annotations_to_qfeatures.R"

rule prepare_annotations:
  input:
   annotated_qf=os.path.join(EXPERIMENT_DIR,"output/output_files/annotated_qf.rds"),
   uniprot_annotated_protein_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/uniprot_annotated_protein_info.txt")
  output:
    conduit_annotations =os.path.join(EXPERIMENT_DIR,"output/output_files/conduit_annotations.tsv")
  log: os.path.join(EXPERIMENT_DIR,"logs/build_conduit/prepare_annotations.log")
  container: "docker://baynec2/conduitr:alpha"
  script: "scripts/prepare_annotations.R"

rule build_conduit:
  input:
    diann_stats= os.path.join(EXPERIMENT_DIR,"output/diann_output/diann.stats.tsv"),
    qfeatures= os.path.join(EXPERIMENT_DIR,"output/output_files/annotated_qf.rds"),
    database=os.path.join(EXPERIMENT_DIR,"input/database_resources/protein_info.txt"),
    annotations= os.path.join(EXPERIMENT_DIR,"output/output_files/conduit_annotations.tsv")
  output:
    conduit = os.path.join(
    EXPERIMENT_DIR,
    "output",
    "output_files",
    f"{config['experiment']}_conduit.rds"
)
  log: os.path.join(EXPERIMENT_DIR,"logs/build_conduit/build_conduit.log")
  container: "docker://baynec2/conduitr:alpha"
  script: "scripts/build_conduit.R"
    
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