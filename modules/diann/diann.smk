import glob
import os
EXPERIMENT_DIR = os.path.join("experiments",config["experiment"])
RAW_FILEPATHS = glob.glob(os.path.join(EXPERIMENT_DIR, "input/raw_files/*.raw"))
#################################################################################
# Generating Spectral Library
#################################################################################
rule generate_diann_spectral_library:
    input:
        fasta = os.path.join(EXPERIMENT_DIR,"input/database_resources/database.fasta"),
        config_file = os.path.join(EXPERIMENT_DIR,"config/generate_diann_spectral_library.cfg")
    output:
        os.path.join(EXPERIMENT_DIR,"input/database_resources/database.predicted.speclib")
    log: os.path.join(EXPERIMENT_DIR,"logs/diann/generate_diann_spectral_library.log")
    container:
        "apptainer/diann2.1.0.sif"
    threads: workflow.cores 
    shell:
        """
        diann --cfg {input.config_file} \
        --fasta {input.fasta} \
        --out-lib experiments/{config[experiment]}/input/database_resources/database \
        --threads {threads} >> {log} 2>&1
        """
################################################################################
# Running DIANN
################################################################################
rule run_diann:
    input:
        raw_files_dir = os.path.join(EXPERIMENT_DIR,"input/raw_files"),
        spectral_library = os.path.join(EXPERIMENT_DIR,"input/database_resources/database.predicted.speclib"),
        fasta = os.path.join(EXPERIMENT_DIR,"input/database_resources/database.fasta"),
        config_file = os.path.join(EXPERIMENT_DIR,"config/run_diann.cfg")
    output:
        report_pr_matrix = os.path.join(EXPERIMENT_DIR,"output/diann_output/diann.pr_matrix.tsv"),
        report_pg_matrix = os.path.join(EXPERIMENT_DIR,"output/diann_output/diann.pg_matrix.tsv")
    log: os.path.join(EXPERIMENT_DIR,"logs/diann/run_diann.log")
    container:
        "apptainer/diann2.1.0.sif"
    threads: workflow.cores 
    shell:
        """
        diann --cfg {input.config_file} \
        --fasta {input.fasta} \
        --out  experiments/{config[experiment]}/output/diann_output/diann \
        --dir {input.raw_files_dir} \
        --lib {input.spectral_library} \
        --threads {threads} --verbose 1 >> {log} 2>&1
        """
################################################################################
# Extracting Detected Proteins
################################################################################
rule extract_detected_proteins:
  input:
    protein_info_df=os.path.join(EXPERIMENT_DIR,"input/database_resources/protein_info.txt"),
    protein_info_fasta =os.path.join(EXPERIMENT_DIR,"input/database_resources/database.fasta"),
    report_pg_matrix=os.path.join(EXPERIMENT_DIR,"output/diann_output/diann.pg_matrix.tsv")
  output:
    detected_protein_info_df = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/detected_protein_info.txt"),
    detected_protein_info_fasta = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/detected_protein.fasta")
  log: os.path.join(EXPERIMENT_DIR,"logs/diann/extract_detected_proteins.log")
  container: "apptainer/conduitR.sif"
  script:
    "scripts/extract_detected_proteins.R"