EXPERIMENT_DIR = os.path.join("experiments",config["experiment"])
################################################################################
# Getting Uniprot Proteome IDs corresponding to the NCBI Taxonomic Identifier
################################################################################        
rule get_uniprot_proteome_ids:
  input:
    os.path.join(EXPERIMENT_DIR,"input/ncbi_taxa_ids.txt")
  output:
    proteome_ids = os.path.join(EXPERIMENT_DIR,"input/proteome_ids.txt")
  log: os.path.join(EXPERIMENT_DIR,"logs/search_space/ncbi_taxonomy/get_uniprot_proteome_ids.log")
  container: "docker://baynec2/conduitr:alpha"
  script:
    "scripts/get_uniprot_proteome_ids.R"