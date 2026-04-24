import glob
import os
EXPERIMENT_DIR = config["experiment_dir"]
RUN_DIR = config["run_dir"]
RAW_FILEPATHS = glob.glob(os.path.join(EXPERIMENT_DIR, "input/raw_files/*.raw"))
#################################################################################
# Generating Spectral Library
#################################################################################
rule generate_diann_spectral_library:
    input:
        fasta = os.path.join(RUN_DIR,"database_resources/database.fasta"),
        config_file = os.path.join(RUN_DIR,"config/diann_spectral_library_base.cfg")
    output:
        os.path.join(RUN_DIR,"database_resources/database.predicted.speclib")
    params:
        out_lib = lambda w, output: os.path.splitext(os.path.splitext(output[0])[0])[0]
    log: os.path.join(RUN_DIR,"logs/diann/generate_diann_spectral_library.log")
    container:
        config["containers"]["diann"]
    threads: workflow.cores
    shell:
        """
        diann --cfg {input.config_file} \
        --fasta {input.fasta} \
        --out-lib {params.out_lib} \
        --met-excision \
        --cut "K*,R*" \
        --missed-cleavages 1 \
        --min-pep-len 7 \
        --max-pep-len 30 \
        --threads {threads} >> {log} 2>&1
        """
################################################################################
# Running DIANN
################################################################################
rule run_diann:
    input:
        raw_files_dir = os.path.join(EXPERIMENT_DIR,"input/raw_files"),
        spectral_library = lambda w: (
            os.path.join(RUN_DIR,"database_resources/database.predicted.speclib")
            if config.get("diann_search_mode", "standard") == "standard"
            else []
        ),
        fasta = os.path.join(RUN_DIR,"database_resources/database.fasta"),
        config_file = os.path.join(RUN_DIR,"config/run_diann.cfg")
    output:
        diann_stats = os.path.join(RUN_DIR,"diann_output/diann.stats.tsv"),
        diann_parquet = os.path.join(RUN_DIR,"diann_output/diann.parquet"),
        diann_pg_matrix = os.path.join(RUN_DIR,"diann_output/diann.pg_matrix.tsv")
    params:
        out = lambda w, output: os.path.join(os.path.dirname(output.diann_stats), "diann"),
        lib_flag = lambda w: (
            f"--lib {os.path.join(RUN_DIR, 'database_resources/database.predicted.speclib')}"
            if config.get("diann_search_mode", "standard") == "standard"
            else ""
        )
    log: os.path.join(RUN_DIR,"logs/diann/run_diann.log")
    container:
        config["containers"]["diann"]
    threads: workflow.cores
    shell:
        """
        diann --cfg {input.config_file} \
        --fasta {input.fasta} \
        --out  {params.out} \
        --dir {input.raw_files_dir} \
        {params.lib_flag} \
        --threads {threads} --verbose 1 >> {log} 2>&1
        """
################################################################################
# Extracting Detected Proteins
################################################################################
rule extract_detected_proteins:
  input:
    protein_info_df=os.path.join(RUN_DIR,"database_resources/protein_info.txt"),
    protein_info_fasta =os.path.join(RUN_DIR,"database_resources/database.fasta"),
    report_pg_matrix=os.path.join(RUN_DIR,"diann_output/diann.pg_matrix.tsv")
  output:
    detected_protein_info_df = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/detected_protein_info.txt"),
    detected_protein_info_fasta = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/detected_protein.fasta")
  log: os.path.join(RUN_DIR,"logs/diann/extract_detected_proteins.log")
  container: config["containers"]["conduitr"]
  script:
    "scripts/extract_detected_proteins.R"
