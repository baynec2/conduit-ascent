EXPERIMENT_DIR = os.path.join("experiments",config["experiment"])

# Get detected protein information from MAGs
rule get_detected_mag_annotations:
  input:
    detected_protein_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/detected_protein_info.txt"),
    mag_annotations = os.path.join(EXPERIMENT_DIR,"input/database_resources/bakta/mag_annotations.txt")
  output:
    bakta_annotated_protein_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/bakta_annotated_protein_info.txt")
  log: os.path.join(EXPERIMENT_DIR,"logs/annotation/MAGs/get_annotations_from_mags.log")
  container: "docker://baynec2/conduitr:alpha"
  script:
    "scripts/get_detected_mag_annotations.R"

# Bakta annotates MAGs via homology to uniref clusters. 
# Our strategy is going to be to get the uniref ids, and use them to annotate the data.
# This way it will be consistent with the rest of the conduit, and traceable.
rule get_annotations_from_uniprot:
  input:
    bakta_annotated_protein_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/bakta_annotated_protein_info.txt")  
  output:
    uniprot_annotated_protein_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/uniprot_annotated_protein_info.txt")
  log: os.path.join(EXPERIMENT_DIR,"logs/annotation/MAGs/get_supplementary_annotations_from_uniprot.log")
  container: "docker://baynec2/conduitr:alpha"
  script: 
    "scripts/get_annotations_from_uniprot.R"

