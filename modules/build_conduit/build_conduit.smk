EXPERIMENT_DIR = config["experiment_dir"]
RUN_DIR = config["run_dir"]
################################################################################
# Processing to QFeatures
################################################################################
rule generate_qfeatures_from_diann_parquet:
  input:
    diann_parquet=os.path.join(RUN_DIR,"diann_output/diann.parquet"),
    sample_annotation=os.path.join(EXPERIMENT_DIR,"input/sample_annotation.txt")
  output:
    qf = os.path.join(RUN_DIR,"output_files/qf.rds")
  log: os.path.join(RUN_DIR,"logs/build_conduit/generate_qfeatures_from_diann_parquet.log")
  container: "docker://baynec2/conduitr:alpha"
  script: "scripts/generate_qfeatures_from_diann_parquet.R"

rule add_annotations_to_qfeatures:
  input:
    qf = os.path.join(RUN_DIR,"output_files/qf.rds"),
    uniprot_annotated_protein_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/uniprot_annotated_protein_info.txt"),
    conduit_annotations = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/conduit_annotations.txt")
  output:
    annotated_qf=os.path.join(RUN_DIR,"output_files/annotated_qf.rds")
  log: os.path.join(RUN_DIR,"logs/build_conduit/add_annotations_to_qfeatures.log")
  container: "docker://baynec2/conduitr:alpha"
  script: "scripts/add_annotations_to_qfeatures.R"

rule build_conduit:
  input:
    diann_stats= os.path.join(RUN_DIR,"diann_output/diann.stats.tsv"),
    qfeatures= os.path.join(RUN_DIR,"output_files/annotated_qf.rds"),
    database=os.path.join(RUN_DIR,"database_resources/protein_info.txt"),
    annotations= os.path.join(RUN_DIR,"database_resources/detected_protein_resources/conduit_annotations.txt"),
    taxonomy = os.path.join(RUN_DIR,"database_resources/taxonomy.txt")
  output:
    conduit = os.path.join(
      RUN_DIR,
      "output_files",
      f"{config['experiment']}_{config['run_name']}_conduit.rds"
    )
  log: os.path.join(RUN_DIR,"logs/build_conduit/build_conduit.log")
  container: "docker://baynec2/conduitr:alpha"
  script: "scripts/build_conduit.R"
