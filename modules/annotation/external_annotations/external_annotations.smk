EXPERIMENT_DIR = os.path.join("experiments",config["experiment"])


# Getting Kegg Information
rule get_kegg_info:
  input:
    uniprot_annotated_protein_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/uniprot_annotated_protein_info.txt")
  output:
    kegg_pathway_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/kegg_pathway_info.txt"),
    kegg_map_pathway_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/kegg_map_pathway_info.txt"),
    kegg_orthology_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/kegg_orthology_info.txt")
  log: os.path.join(EXPERIMENT_DIR,"logs/annotation/get_kegg_info.log")
  container: "docker://baynec2/conduitr:alpha"
  script:
    "scripts/get_kegg_info.R"

# Getting Pfam info
rule get_pfam_resources:
    output:
        "resources/annotation/pfam/pfamA.txt"
    log:
        "resources/annotation/pfam/get_pfam_resources.log"
    container:
        "docker://baynec2/conduitr:alpha"
    shell:
        """
        mkdir -p resources/annotation/pfam
        curl -L -o resources/annotation/pfam/pfamA.txt.gz \
            https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/pfamA.txt.gz \
            &> {log}
        gunzip -c resources/annotation/pfam/pfamA.txt.gz > {output} 2>> {log}
        """

# Getting Pfam info
rule get_pfam_info:
  input:
    uniprot_annotated_protein_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/uniprot_annotated_protein_info.txt"),
    pfam_db ="resources/annotation/pfam/pfamA.txt"
  output:
    pfam_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/pfam_info.txt")
  log:
    os.path.join(EXPERIMENT_DIR,"logs/annotation/ncbi_taxonomy/pfam_info.log")
  container:
    "docker://baynec2/conduitr:alpha"
  script:  
    "scripts/get_pfam_info.R"

rule get_cazy_resource:
    output:
        "resources/annotation/cazy/CAZyDB.07302020.fam-activities.txt"
    log:
        os.path.join(EXPERIMENT_DIR, "logs/annotation/cazy/get_cazy_resource.log")
    container:
        "docker://baynec2/conduitr:alpha"
    shell:
        """
        mkdir -p resources/annotation/cazy
        curl -L -o {output} https://bcb.unl.edu/dbCAN2/download/Databases/CAZyDB.07302020.fam-activities.txt \
            &> {log}
        """
# Getting Cazyme Information
rule get_cazy_info:
  input:
    uniprot_annotated_protein_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/uniprot_annotated_protein_info.txt"),
    cazy_resource = "resources/annotation/cazy/CAZyDB.07302020.fam-activities.txt"
  output:
    cazy_class_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/cazy_class_info.txt"),
    cazy_family_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/cazy_family_info.txt")
  log: os.path.join(EXPERIMENT_DIR,"logs/annotation/get_cazy_info.log")
  container: "docker://baynec2/conduitr:alpha"
  script:
    "scripts/get_cazy_info.R"

rule get_eggnog_resources:
    output:
        eggnog_resource = "resources/annotation/eggnog/e5.og_annotations.tsv"
    log: "resources/annotation/eggnog/get_eggnog_resources.log"
    container: "docker://baynec2/conduitr:alpha"
    shell:
        """
        mkdir -p resources/annotation/eggnog
        curl -L -o {output.eggnog_resource} http://eggnog6.embl.de/download/eggnog_5.0/e5.og_annotations.tsv \
            &> {log}
        """
# Getting Eggnog Information
rule get_eggnog_info:
  input:
    uniprot_annotated_protein_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/uniprot_annotated_protein_info.txt"),
    eggnog_resource = "resources/annotation/eggnog/e5.og_annotations.tsv"
  output:
    eggnog_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/eggnog_info.txt"),
    eggnog_code_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/eggnog_code_info.txt")
  log: os.path.join(EXPERIMENT_DIR,"logs/annotation/get_eggnog_info.log")
  container: "docker://baynec2/conduitr:alpha"
  script:
    "scripts/get_eggnog_info.R"

# Getting GO information (from UniProt annotations)
rule get_go_info:
  input:
    uniprot_annotated_protein_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/uniprot_annotated_protein_info.txt"),
  output:
    go_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/go_info.txt")
  log: os.path.join(EXPERIMENT_DIR,"logs/annotation/get_go_info.log")
  container: "docker://baynec2/conduitr:alpha"
  script:
    "scripts/get_go_info.R"

# Consolidating annotations
rule consolidate_annotations:
  input:
    # Annotations
    go_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/go_info.txt"),
    kegg_pathway_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/kegg_pathway_info.txt"),
    kegg_map_pathway_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/kegg_map_pathway_info.txt"),
    kegg_orthology_info =  os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/kegg_orthology_info.txt"),
    pfam_info =  os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/pfam_info.txt"),
    cazy_class_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/cazy_class_info.txt"),
    cazy_family_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/cazy_family_info.txt"),
    eggnog_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/eggnog_info.txt"),
    eggnog_code_info = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/eggnog_code_info.txt"),
    # QF
    qf = os.path.join(EXPERIMENT_DIR,"output/output_files/qf.rds")
  output:
    conduit_annotations = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/conduit_annotations.txt")
  log: os.path.join(EXPERIMENT_DIR,"logs/annotation/consolidate_annotations.log")
  container: "docker://baynec2/conduitr:alpha"
  script:
    "scripts/consolidate_annotations.R"
