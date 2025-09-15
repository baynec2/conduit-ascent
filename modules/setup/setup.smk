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