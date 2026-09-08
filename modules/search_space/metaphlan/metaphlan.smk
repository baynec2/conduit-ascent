import glob
import os
EXPERIMENT_DIR = config["experiment_dir"]
RUN_DIR = config["run_dir"]
METAPHLAN_DB_DIR = config["metaphlan_database_dir"]

################################################################################
# MetaPhlAn Database Management
################################################################################
# `mpa_latest` is the canonical sentinel file emitted by `metaphlan --install`.
# Using a real file as the rule output (not the whole directory) lets us point
# metaphlan_database_dir at an existing off-repo DB on a machine profile
# (no symlinks; apptainer binds the parent path). When the file is present
# snakemake skips the install; when it isn't, the install rule populates it.
rule download_metaphlan_resources:
    output:
        marker_index = os.path.join(METAPHLAN_DB_DIR, "mpa_latest")
    container: config["containers"]["metaphlan"]
    log: os.path.join(METAPHLAN_DB_DIR, "logs/download_metaphlan_resources.log")
    shell: "metaphlan --install --db_dir " + METAPHLAN_DB_DIR + " 2> {log}"

###############################################################################
# Generating MetaPhlAn Output Files
###############################################################################
rule run_metaphlan:
    input:
        fastq=os.path.join(EXPERIMENT_DIR, "input/fastq_files/{sample}.fastq.gz"),
        marker_index=os.path.join(METAPHLAN_DB_DIR, "mpa_latest")
    output:
        profile=os.path.join(RUN_DIR, "metaphlan/{sample}_profile.txt"),
        mapout=os.path.join(RUN_DIR, "metaphlan/{sample}.mapout.txt")
    container: config["containers"]["metaphlan"]
    benchmark:
        os.path.join(RUN_DIR, "benchmarks/search_space/metaphlan/run_metaphlan_{sample}.tsv")
    log:
        os.path.join(RUN_DIR, "logs/search_space/metaphlan/run_metaphlan_{sample}.log")
    threads: workflow.cores
    shell:
        """
        metaphlan {input.fastq} \
            --input_type fastq \
            --nproc {threads} \
            --db_dir """ + METAPHLAN_DB_DIR + """ \
            --mapout {output.mapout} \
            -o {output.profile} \
            >> {log} 2>&1
        """
###############################################################################
# Merging MetaPhlAn Profile Files
###############################################################################
rule merge_profiles:
    input:
        metaphlan_profiles = expand(os.path.join(RUN_DIR, "metaphlan/{sample}_profile.txt"),
                         sample=glob_wildcards(os.path.join(EXPERIMENT_DIR, "input/fastq_files/{sample}.fastq.gz")).sample)
    output:
        merged_profiles = os.path.join(RUN_DIR, "metaphlan/merged_profiles.txt")
    container: config["containers"]["conduitr"]
    log: os.path.join(RUN_DIR, "logs/search_space/metaphlan/combine_metaphlan_output.log")
    script: "scripts/merge_profiles.R"
###############################################################################
# Applying Threshold and Converting to ncbi_taxonomy_id input
###############################################################################
rule call_ncbi_taxa_ids:
    input:
        merged_profiles = os.path.join(RUN_DIR, "metaphlan/merged_profiles.txt")
    output:
        ncbi_taxa_ids = os.path.join(RUN_DIR, "ncbi_taxa_ids.txt")
    container: config["containers"]["conduitr"]
    log: os.path.join(RUN_DIR, "logs/search_space/metaphlan/call_ncbi_taxa_ids.log")
    script:
        "scripts/call_ncbi_taxa_ids.R"