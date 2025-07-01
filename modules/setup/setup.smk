import os
import logging

# Define a function to setup logging for a rule
def setup_rule_logging(log_file):
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='[%(asctime)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger()

EXPERIMENT_DIR = os.path.join("experiments",config["experiment"])
################################################################################
# Configuration Setup Rules
################################################################################
# Handle DIANN spectral library config
rule setup_diann_spectral_library_config:
    output:
        output_config_file = os.path.join(EXPERIMENT_DIR,"config/generate_diann_spectral_library.cfg")
    params:
        default_config = "config/generate_diann_spectral_library.cfg",
        selected_config = config.get("generate_diann_spectral_library_config")
    log: os.path.join(EXPERIMENT_DIR,"logs/setup/setup_diann_spectral_library_config.log")
    run:
        logger = setup_rule_logging(log[0])
        logger.info("Starting setup_diann_spectral_library_config rule")

        if params.selected_config == params.default_config:
            shell("cp {params.default_config} {output.output_config_file}")
            logger.info("Copied default spectral library config to experiment directory")
        else:
            shell("cp {params.selected_config} {output.output_config_file}")
            logger.info("Copied custom spectral library config to experiment directory")

        logger.info("setup_diann_spectral_library_config rule completed")

# Handle DIANN run config
rule setup_diann_run_config:
    output:
        diann_run_config_file = os.path.join(EXPERIMENT_DIR,"config/run_diann.cfg")
    params:
        default_config = "config/run_diann.cfg",
        selected_config = config.get("run_diann_config")
    log: os.path.join(EXPERIMENT_DIR,"logs/setup/setup_diann_run_config.log")
    run:
        logger = setup_rule_logging(log[0])
        logger.info("Starting setup_diann_run_config rule")
                
        if params.selected_config == params.default_config:
            shell("cp {params.default_config} {output.diann_run_config_file}")
            logger.info("Copied default run config to experiment directory")
        else:
            shell("cp {params.selected_config} {output.diann_run_config_file}")
            logger.info("Copied custom run config to experiment directory")
            
        logger.info("setup_diann_run_config rule completed")

################################################################################
# Apptainer Setup Rules
################################################################################
# Checking if apptainer is installed, if using apptainer
rule check_apptainer:
    output:
         os.path.join(EXPERIMENT_DIR,"logs/setup/apptainer_checked")
    log: os.path.join(EXPERIMENT_DIR,"logs/setup/check_apptainer.log")
    run:
        logger = setup_rule_logging(log[0])
        logger.info("Starting Apptainer check...")
        try:
            shell("apptainer --version")
            logger.info(f"Apptainer was found sucessfully")
        except Exception as e:
            logger.warning("Apptainer is not installed or not found in PATH.")
            logger.warning("Please install Apptainer before running this workflow if apptainer is desired.")
            logger.warning("Installation instructions: https://apptainer.org/docs/admin/main/installation.html")
        shell("touch {output}")
        logger.info("check_apptainer rule completed successfully.")

# Build conduitR apptainer 
rule build_conduitR_apptainer:
    input:
        "apptainer/conduitR.def"
    output:
        "apptainer/conduitR.sif"
    log: os.path.join(EXPERIMENT_DIR,"logs/setup/build_conduitR_apptainer.log")
    run:
        logger = setup_rule_logging(log[0])
        logger.info("Starting build_conduitR_apptainer rule")
        try:
            logger.info("Attempting to build conduitR apptainer...")
            shell("apptainer build {output} {input}")
            logger.info("conduitR apptainer successfully built")
        except Exception as e:
            logger.warning(f"Unable to build conduitR apptainer: {str(e)}")
            logger.warning("If apptainer is desired, please ensure it is properly installed.")
            # Add a blank .sif file to the output directory to allow the workflow to continue to run
            shell("touch {output}")
        logger.info("build_conduitR_apptainer rule completed")

# Build diann apptainer 
rule build_diann_apptainer:
    input:
        "apptainer/diann2.1.0.def"
    output:
        "apptainer/diann2.1.0.sif"
    log: os.path.join(EXPERIMENT_DIR,"logs/setup/build_diann_apptainer.log")
    run:
        logger = setup_rule_logging(log[0])
        try:
            logger.info("Attempting to build diann2.1.0 apptainer...")
            shell("apptainer build {output} {input}")
            logger.info("diann2.1.0 apptainer successfully built")
        except Exception as e:
            logger.warning(f"Unable to build diann2.1.0 apptainer: {str(e)}")
            logger.warning("If apptainer is desired, please ensure it is properly installed.")
            # Add a blank .sif file to the output directory to allow the workflow to continue to run
            shell("touch {output}")
        logger.info("build_diann_apptainer rule completed") 