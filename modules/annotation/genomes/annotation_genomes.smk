EXPERIMENT_DIR = config["experiment_dir"]
RUN_DIR = config["run_dir"]

# Get detected protein information from genome annotations
rule get_detected_genome_annotations:
  input:
    detected_protein_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/detected_protein_info.txt"),
    genome_annotations = os.path.join(RUN_DIR,"database_resources/genome_annotations.txt")
  output:
    bakta_annotated_protein_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/bakta_annotated_protein_info.txt")
  log: os.path.join(RUN_DIR,"logs/annotation/genomes/get_annotations_from_genomes.log")
  container: config["containers"]["conduitr"]
  script:
    "scripts/get_detected_genome_annotations.R"

# Bakta annotates genomes via homology to uniref clusters.
# Our strategy is going to be to get the uniref ids, and use them to annotate the data.
# This way it will be consistent with the rest of the conduit, and traceable.
rule get_annotations_from_uniprot:
  input:
    bakta_annotated_protein_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/bakta_annotated_protein_info.txt")
  output:
    uniprot_annotated_protein_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/uniprot_annotated_protein_info.txt")
  log: os.path.join(RUN_DIR,"logs/annotation/genomes/get_supplementary_annotations_from_uniprot.log")
  # UniProt-API annotation is network-bound; cap workers (and thus concurrent
  # API requests) to a modest count. snakemake@threads is propagated into the
  # script to size the parallel worker pool (see the script header).
  threads: min(8, workflow.cores)
  container: config["containers"]["conduitr"]
  script:
    "scripts/get_annotations_from_uniprot.R"
