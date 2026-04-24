import os

EXPERIMENT_DIR = config["experiment_dir"]
RUN_DIR = config["run_dir"]

localrules: setup_diann_spectral_library_config, setup_diann_run_config, write_run_manifest

################################################################################
# Configuration Setup Rules
################################################################################
# Handle DIANN spectral library config
rule setup_diann_spectral_library_config:
    output:
        output_config_file = os.path.join(RUN_DIR,"config/generate_diann_spectral_library.cfg")
    params:
        selected_config = config.get("generate_diann_spectral_library_config")
    log: os.path.join(RUN_DIR,"logs/setup/setup_diann_spectral_library_config.log")
    run:
        import shutil
        shutil.copy(params.selected_config, output.output_config_file)

# Handle DIANN run config
rule setup_diann_run_config:
    output:
        diann_run_config_file = os.path.join(RUN_DIR,"config/run_diann.cfg")
    params:
        selected_config = (
            "config/run_diann_infinidia.cfg"
            if config.get("diann_search_mode") == "infinidia"
            else config.get("run_diann_config")
        )
    log: os.path.join(RUN_DIR,"logs/setup/setup_diann_run_config.log")
    run:
        import shutil
        shutil.copy(params.selected_config, output.diann_run_config_file)

# Point-in-time reproducibility snapshot: resolved merged config + git SHA +
# sample_annotation hash + cfg-file hashes + container tags. See
# modules/setup/scripts/write_run_manifest.py for the recorded fields.
rule write_run_manifest:
    input:
        sample_annotation = os.path.join(EXPERIMENT_DIR, config["sample_annotation"])
    output:
        manifest = os.path.join(RUN_DIR, "manifest.json")
    params:
        repo_root = workflow.basedir
    log: os.path.join(RUN_DIR, "logs/setup/write_run_manifest.log")
    script: "scripts/write_run_manifest.py"