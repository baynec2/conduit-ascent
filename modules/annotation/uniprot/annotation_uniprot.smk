EXPERIMENT_DIR = config["experiment_dir"]
RUN_DIR = config["run_dir"]

# Get detected protein information from Uniprot
rule get_annotations_from_uniprot:
  input:
    detected_protein_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/detected_protein_info.txt")
  output:
    uniprot_annotated_protein_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/uniprot_annotated_protein_info.txt")
  benchmark: os.path.join(RUN_DIR,"benchmarks/annotation/uniprot/get_annotations_from_uniprot.tsv")
  log: os.path.join(RUN_DIR,"logs/annotation/uniprot/get_annotations_from_uniprot.log")
  # conduitR::get_annotations_from_uniprot() parallelises UniProt API batches;
  # the script caps its worker pool to this thread count (see the R script).
  threads: min(8, workflow.cores)
  container: config["containers"]["conduitr"]
  script:
    "scripts/get_annotations_from_uniprot.R"
