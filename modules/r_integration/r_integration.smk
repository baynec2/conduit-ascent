EXPERIMENT_DIR = os.path.join("experiments",config["experiment"])
################################################################################
# Ingestion into R
################################################################################
# Prepare a QFeatures object with all the data
rule prepare_qf:
  input:
    annotation = os.path.join(EXPERIMENT_DIR,"input/sample_annotation.txt"),
    precursor_matrix=os.path.join(EXPERIMENT_DIR,"output/output_files/precursor_matrix.tsv"),
    peptide_matrix = os.path.join(EXPERIMENT_DIR,"output/output_files/peptide_matrix.tsv"),
    protein_group_matrix=os.path.join(EXPERIMENT_DIR,"output/output_files/protein_group_matrix.tsv"),
    domain_matrix=os.path.join(EXPERIMENT_DIR,"output/output_files/domain_matrix.tsv"),
    kingdom_matrix=os.path.join(EXPERIMENT_DIR,"output/output_files/kingdom_matrix.tsv"),
    phylum_matrix=os.path.join(EXPERIMENT_DIR,"output/output_files/phylum_matrix.tsv"),
    class_matrix=os.path.join(EXPERIMENT_DIR,"output/output_files/class_matrix.tsv"),
    order_matrix=os.path.join(EXPERIMENT_DIR,"output/output_files/order_matrix.tsv"),
    family_matrix=os.path.join(EXPERIMENT_DIR,"output/output_files/family_matrix.tsv"),
    genus_matrix=os.path.join(EXPERIMENT_DIR,"output/output_files/genus_matrix.tsv"),
    species_matrix=os.path.join(EXPERIMENT_DIR,"output/output_files/species_matrix.tsv"),
    go_matrix = os.path.join(EXPERIMENT_DIR,"output/output_files/go_matrix.tsv"),
    go_taxa_matrix = os.path.join(EXPERIMENT_DIR,"output/output_files/go_taxa_matrix.tsv"),
    subcellular_locations_matrix = os.path.join(EXPERIMENT_DIR,"output/output_files/subcellular_locations_matrix.tsv"),
    kegg_pathway_matrix = os.path.join(EXPERIMENT_DIR,"output/output_files/kegg_pathway_matrix.tsv"),
    kegg_ko_matrix = os.path.join(EXPERIMENT_DIR,"output/output_files/kegg_ko_matrix.tsv")
  output:
    qf= os.path.join(EXPERIMENT_DIR,"output/output_files/qf.rds")
  log: os.path.join(EXPERIMENT_DIR,"logs/r_integration/prepare_qf.log")
  container: "docker://baynec2/conduitr:alpha"
  script:
    "scripts/prepare_qf.R"

# prepare metrics    
rule extract_metrics:
  input:
    protein_info = os.path.join(EXPERIMENT_DIR,"output/database_resources/protein_info.txt"),
    detected_protein_info = os.path.join(EXPERIMENT_DIR,"output/database_resources/detected_protein_resources/detected_protein_info.txt")
  output:
    database_taxonomy=os.path.join(EXPERIMENT_DIR,"output/output_files/database_taxonomy.tsv"),
    database_metrics=os.path.join(EXPERIMENT_DIR,"output/output_files/database_metrics.tsv"),
    detected_protein_taxonomy = os.path.join(EXPERIMENT_DIR,"output/output_files/detected_protein_taxonomy.tsv"),
    detected_protein_metrics = os.path.join(EXPERIMENT_DIR,"output/output_files/detected_protein_metrics.tsv"),
    combined_metrics = os.path.join(EXPERIMENT_DIR,"output/output_files/combined_metrics.tsv")
  log: os.path.join(EXPERIMENT_DIR,"logs/r_integration/extract_metrics.log")
  container: "docker://baynec2/conduitr:alpha"
  script:
    "scripts/extract_metrics.R"
    
# Prepare a conduit object
# Experiment metadata = date of experiment
# Database taxonomy = n proteins belonging to each taxa in database
# database metrics = n_proteins, n_organism_types, domain, kingdom, phylum, class, order, family, genus,species
# detected protein taxonomy = n proteins detected belonging to each taxa
# detected protein metrics = n_protein_groups,n_uniquely_ided_proteins,n_precursors,n_peptides n_organism_types, domain, kingdom, phylum, class, order, family, genus,species
rule create_conduit:
  input:
    qf = os.path.join(EXPERIMENT_DIR,"output/output_files/qf.rds"),
    database_taxonomy = os.path.join(EXPERIMENT_DIR,"output/output_files/database_taxonomy.tsv"),
    combined_metrics = os.path.join(EXPERIMENT_DIR,"output/output_files/combined_metrics.tsv"),
    detected_protein_taxonomy = os.path.join(EXPERIMENT_DIR,"output/output_files/detected_protein_taxonomy.tsv")
  output:
    conduit_obj = os.path.join(EXPERIMENT_DIR,"output/output_files/conduit_output.rds")
  log: os.path.join(EXPERIMENT_DIR,"logs/r_integration/create_conduit.log")
  container: "docker://baynec2/conduitr:alpha"
    # From this point- the conduit object contains all the stuff you need. 
    # Load that into Conduit GUI.
  script:
    "scripts/create_conduit_obj.R"