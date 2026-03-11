EXPERIMENT_DIR = os.path.join("experiments",config["experiment"])

# Get detected protein information from Uniprot
rule get_annotations_from_uniprot:
  input:
    detected_protein_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/detected_protein_info.txt")
  output:
    uniprot_annotated_protein_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/uniprot_annotated_protein_info.txt")
  log: os.path.join(EXPERIMENT_DIR,"logs/annotation/uniprot/get_annotations_from_uniprot.log")
  container: "docker://baynec2/conduitr:alpha"
  script:
    "scripts/get_annotations_from_uniprot.R"
