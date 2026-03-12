EXPERIMENT_DIR = config["experiment_dir"]
RUN_DIR = config["run_dir"]
################################################################################
# Getting Uniprot Proteome IDs corresponding to the NCBI Taxonomic Identifier
################################################################################
# For ncbi_taxonomy_id and peptidotyping methods ncbi_taxa_ids.txt is user-provided
# in EXPERIMENT_DIR/input/. A seed rule copies it into RUN_DIR so that all
# downstream rules can always read from RUN_DIR consistently.
# For the metaphlan method ncbi_taxa_ids.txt is generated directly into RUN_DIR
# by the metaphlan module.
if config["search_space_method"] in ("ncbi_taxonomy_id", "peptidotyping"):
    rule seed_ncbi_taxa_ids:
        input: os.path.join(EXPERIMENT_DIR, "input/ncbi_taxa_ids.txt")
        output: os.path.join(RUN_DIR, "ncbi_taxa_ids.txt")
        log: os.path.join(RUN_DIR, "logs/search_space/ncbi_taxonomy/seed_ncbi_taxa_ids.log")
        shell: "cp {input} {output} 2> {log}"

rule get_uniprot_proteome_ids:
  input:
    os.path.join(RUN_DIR, "ncbi_taxa_ids.txt")
  output:
    proteome_ids = os.path.join(RUN_DIR,"proteome_ids.txt")
  log: os.path.join(RUN_DIR,"logs/search_space/ncbi_taxonomy/get_uniprot_proteome_ids.log")
  container: "docker://baynec2/conduitr:alpha"
  script:
    "scripts/get_uniprot_proteome_ids.R"
