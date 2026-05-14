EXPERIMENT_DIR = config["experiment_dir"]
RUN_DIR = config["run_dir"]


def search_space_detection_inputs(wildcards):
    """Optional taxon-detection artifacts for build_conduit, keyed by
    metric-slot name. Methods that don't produce a detection artifact return
    an empty dict, so no extra inputs are added to the rule's DAG."""
    m = config["search_space_method"]
    base = os.path.join(RUN_DIR, "database_resources")
    if m == "unipept_peptidotyping":
        return {
            "peptidotyping_first_pass":  os.path.join(base, "peptidotyping/first_pass_fdr_results.tsv"),
            "peptidotyping_second_pass": os.path.join(base, "peptidotyping/second_pass_fdr_results.tsv"),
        }
    if m == "genome_peptidotyping":
        return {
            "peptidotyping_first_pass":  os.path.join(base, "genome_peptidotyping/first_pass_fdr_results.tsv"),
            "peptidotyping_second_pass": os.path.join(base, "genome_peptidotyping/second_pass_fdr_results.tsv"),
        }
    if m == "unipept_hapid":
        return {"hapid_greedy_selection": os.path.join(base, "unipept_hapid/unipept_hapid_greedy_selection.tsv")}
    if m == "hapid":
        return {"hapid_greedy_selection": os.path.join(base, "hapid/hapid_greedy_selection.tsv")}
    return {}


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
  container: config["containers"]["conduitr"]
  script: "scripts/generate_qfeatures_from_diann_parquet.R"

rule add_annotations_to_qfeatures:
  input:
    qf = os.path.join(RUN_DIR,"output_files/qf.rds"),
    uniprot_annotated_protein_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/uniprot_annotated_protein_info.txt"),
    conduit_annotations = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/conduit_annotations.txt"),
    protein_info = os.path.join(RUN_DIR,"database_resources/protein_info.txt")
  output:
    annotated_qf=os.path.join(RUN_DIR,"output_files/annotated_qf.rds")
  log: os.path.join(RUN_DIR,"logs/build_conduit/add_annotations_to_qfeatures.log")
  container: config["containers"]["conduitr"]
  script: "scripts/add_annotations_to_qfeatures.R"

rule build_conduit:
  input:
    unpack(search_space_detection_inputs),
    diann_stats= os.path.join(RUN_DIR,"diann_output/diann.stats.tsv"),
    qfeatures= os.path.join(RUN_DIR,"output_files/annotated_qf.rds"),
    database=os.path.join(RUN_DIR,"database_resources/protein_info.txt"),
    annotations= os.path.join(RUN_DIR,"database_resources/detected_protein_resources/conduit_annotations.txt"),
    taxonomy = os.path.join(RUN_DIR,"database_resources/taxonomy.txt"),
    diann_spectral_lib_config = config["diann_spectral_library_base_config"],
    diann_run_config          = config["run_diann_config"]
  params:
    workflow_version  = open("VERSION").read().strip(),
    snakemake_version = __import__('snakemake').__version__
  output:
    conduit = os.path.join(
      RUN_DIR,
      "output_files",
      f"{config['experiment']}_{config['run_name']}_conduit.rds"
    )
  log: os.path.join(RUN_DIR,"logs/build_conduit/build_conduit.log")
  container: config["containers"]["conduitr"]
  script: "scripts/build_conduit.R"
