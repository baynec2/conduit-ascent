EXPERIMENT_DIR = os.path.join("experiments",config["experiment"])
# Get detected protein information from Uniprot
rule get_annotations_from_uniprot:
  input:
    detected_protein_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/detected_protein_info.txt")
  output:
    uniprot_annotated_protein_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/uniprot_annotated_protein_info.txt")
  log: os.path.join(EXPERIMENT_DIR,"logs/annotation/ncbi_taxonomy/get_annotations_from_uniprot.log")
  container: "apptainer/conduitR.sif"
  script:
    "scripts/get_annotations_from_uniprot.R"
    
# Extracting GO infromation from data frame
rule extract_go_info:
    input:
     uniprot_annotated_protein_info=os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/uniprot_annotated_protein_info.txt")
    output:
     go_annotations=os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/go_annotations.txt")
    log: os.path.join(EXPERIMENT_DIR,"logs/annotation/ncbi_taxonomy/extract_go_info.log")
    container: "apptainer/conduitR.sif"
    script:
     "scripts/extract_go_info.R" 

# Extracting cellular location information
rule extract_cellular_location_info:
    input:
      uniprot_annotated_protein_info=os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/uniprot_annotated_protein_info.txt")
    output:
      subcellular_locations=os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/subcellular_locations.txt")
    log: os.path.join(EXPERIMENT_DIR,"logs/annotation/ncbi_taxonomy/extract_cellular_location_info.log")
    container: "apptainer/conduitR.sif"
    script:
     "scripts/extract_subcellular_locations.R" 
     
# Getting all Kegg infromation 
rule get_kegg_info:
    input:
      uniprot_annotated_protein_info=os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/uniprot_annotated_protein_info.txt")
    output:
       kegg_annotations=os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/kegg_annotations.txt")
    log: os.path.join(EXPERIMENT_DIR,"logs/annotation/ncbi_taxonomy/get_kegg_info.log")
    container: "apptainer/conduitR.sif"
    script:
     "scripts/get_kegg_info.R" 
     
# Getting all taxonomy infromation 
rule extract_detected_taxonomy:
    input:
      uniprot_annotated_protein_info=os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/uniprot_annotated_protein_info.txt")
    output:
       detected_taxonomy=os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/detected_taxonomy.txt")
    log: os.path.join(EXPERIMENT_DIR,"logs/annotation/ncbi_taxonomy/extract_detected_taxonomy.log")
    container: "apptainer/conduitR.sif"
    script:
     "scripts/extract_detected_taxonomy.R"      
  