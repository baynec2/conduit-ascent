import glob
import os
EXPERIMENT_DIR = config["experiment_dir"]
RUN_DIR = config["run_dir"]

################################################################################
# MetaPhlAn Database Management
################################################################################
rule download_metaphlan_resources:
    output:
        database_dir = directory("resources/metaphlan")
    container: config["containers"]["metaphlan"]
    log: "resources/metaphlan/logs/download_metaphlan_resources.log"
    shell: "metaphlan --install --db_dir {output.database_dir} 2> {log}"

###############################################################################
# Generating MetaPhlAn Output Files
###############################################################################
rule run_metaphlan:
    input:
        fastq=os.path.join(EXPERIMENT_DIR, "input/fastq_files/{sample}.fastq.gz")
    output:
        profile=os.path.join(RUN_DIR, "metaphlan/{sample}_profile.txt")
    container: config["containers"]["metaphlan"]
    log:
        os.path.join(RUN_DIR, "logs/search_space/metaphlan/run_metaphlan_{sample}.log")
    threads: workflow.cores
    shell:
        """
        metaphlan {input.fastq} \
            --input_type fastq \
            --nproc {threads} \
            --db_dir resources/metaphlan/ \
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