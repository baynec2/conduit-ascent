EXPERIMENT_DIR = config["experiment_dir"]
RUN_DIR = config["run_dir"]

PFAM_DB_URL = config.get(
    "pfam_db_url",
    "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/pfamA.txt.gz",
)
CAZY_DB_URL = config.get(
    "cazy_db_url",
    "https://pro.unl.edu/dbCAN2/download/Databases/CAZyDB.07302020.fam-activities.txt",
)
EGGNOG_DB_URL = config.get(
    "eggnog_db_url",
    "http://eggnog6.embl.de/download/eggnog_5.0/e5.og_annotations.tsv",
)

# Authoritative term-name dictionaries used to fill `description` for the
# eggNOG-mapper-derived annotation types (go/kegg_*/pfam/ec_number/brite).
# Pinned to specific releases for reproducibility; bump intentionally.
DICT_DIR = "resources/annotation/dictionaries"
GO_OBO_URL = config.get(
    "go_obo_url",
    "http://release.geneontology.org/2026-05-19/ontology/go-basic.obo",
)
KEGG_REST_URL = config.get("kegg_rest_url", "https://rest.kegg.jp")
ENZYME_DAT_URL = config.get(
    "enzyme_dat_url",
    "https://ftp.expasy.org/databases/enzyme/enzyme.dat",
)
PFAM_CLANS_URL = config.get(
    "pfam_clans_url",
    "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/Pfam-A.clans.tsv.gz",
)

# dbCAN's older bcb.unl.edu host has an incomplete cert chain in some containers.
CAZY_DB_URL = CAZY_DB_URL.replace(
    "https://bcb.unl.edu/dbCAN2/",
    "https://pro.unl.edu/dbCAN2/",
).replace(
    "http://bcb.unl.edu/dbCAN2/",
    "https://pro.unl.edu/dbCAN2/",
)


# Getting Kegg Information
rule get_kegg_info:
  input:
    uniprot_annotated_protein_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/uniprot_annotated_protein_info.txt")
  output:
    kegg_pathway_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/kegg_pathway_info.txt"),
    kegg_map_pathway_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/kegg_map_pathway_info.txt"),
    kegg_orthology_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/kegg_orthology_info.txt")
  log: os.path.join(RUN_DIR,"logs/annotation/get_kegg_info.log")
  container: config["containers"]["conduitr"]
  script:
    "scripts/get_kegg_info.R"

# Getting Pfam info
rule get_pfam_resources:
    output:
        "resources/annotation/pfam/pfamA.txt"
    log:
        "resources/annotation/pfam/get_pfam_resources.log"
    container:
        config["containers"]["conduitr"]
    params:
        pfam_db_url = PFAM_DB_URL
    shell:
        """
        mkdir -p resources/annotation/pfam
        curl -L -o resources/annotation/pfam/pfamA.txt.gz \
            {params.pfam_db_url} \
            &> {log}
        gunzip -c resources/annotation/pfam/pfamA.txt.gz > {output} 2>> {log}
        """

# Getting Pfam info
rule get_pfam_info:
  input:
    uniprot_annotated_protein_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/uniprot_annotated_protein_info.txt"),
    pfam_db ="resources/annotation/pfam/pfamA.txt"
  output:
    pfam_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/pfam_info.txt")
  log:
    os.path.join(RUN_DIR,"logs/annotation/ncbi_taxonomy/pfam_info.log")
  container:
    config["containers"]["conduitr"]
  script:
    "scripts/get_pfam_info.R"

rule get_cazy_resource:
    output:
        "resources/annotation/cazy/cazy_db.txt"
    log:
        os.path.join(RUN_DIR, "logs/annotation/cazy/get_cazy_resource.log")
    container:
        config["containers"]["conduitr"]
    params:
        cazy_db_url = CAZY_DB_URL
    shell:
        """
        mkdir -p resources/annotation/cazy
        if ! curl --fail --location -o {output} {params.cazy_db_url} &> {log}; then
            printf '%s\n' 'Retrying CAZy download without certificate verification.' >> {log}
            rm -f {output}
            curl --fail --location --insecure -o {output} {params.cazy_db_url} \
                &>> {log}
        fi
        """
# Getting Cazyme Information
rule get_cazy_info:
  input:
    uniprot_annotated_protein_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/uniprot_annotated_protein_info.txt"),
    cazy_resource = "resources/annotation/cazy/cazy_db.txt"
  output:
    cazy_class_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/cazy_class_info.txt"),
    cazy_family_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/cazy_family_info.txt")
  log: os.path.join(RUN_DIR,"logs/annotation/get_cazy_info.log")
  container: config["containers"]["conduitr"]
  script:
    "scripts/get_cazy_info.R"

rule get_eggnog_resources:
    output:
        eggnog_resource = "resources/annotation/eggnog/e5.og_annotations.tsv"
    log: "resources/annotation/eggnog/get_eggnog_resources.log"
    container: config["containers"]["conduitr"]
    params:
        eggnog_db_url = EGGNOG_DB_URL
    shell:
        """
        mkdir -p resources/annotation/eggnog
        curl -L -o {output.eggnog_resource} {params.eggnog_db_url} \
            &> {log}
        """
# Getting Eggnog Information
rule get_eggnog_info:
  input:
    uniprot_annotated_protein_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/uniprot_annotated_protein_info.txt"),
    eggnog_resource = "resources/annotation/eggnog/e5.og_annotations.tsv"
  output:
    eggnog_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/eggnog_info.txt"),
    eggnog_code_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/eggnog_code_info.txt")
  log: os.path.join(RUN_DIR,"logs/annotation/get_eggnog_info.log")
  container: config["containers"]["conduitr"]
  script:
    "scripts/get_eggnog_info.R"

# Getting GO information (from UniProt annotations)
rule get_go_info:
  input:
    uniprot_annotated_protein_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/uniprot_annotated_protein_info.txt"),
  output:
    go_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/go_info.txt")
  log: os.path.join(RUN_DIR,"logs/annotation/get_go_info.log")
  container: config["containers"]["conduitr"]
  script:
    "scripts/get_go_info.R"

# Downloading authoritative term-name dictionaries (once, cached in resources/).
# These name the eggNOG-mapper-derived accessions (GO/KEGG/EC/Pfam/BRITE) so the
# `description` column is populated from an explicit external source rather than
# borrowed from the UniProt-side annotations.
rule download_annotation_dictionaries:
    output:
        go_obo       = os.path.join(DICT_DIR, "go-basic.obo"),
        kegg_ko      = os.path.join(DICT_DIR, "kegg_ko.tsv"),
        kegg_pathway = os.path.join(DICT_DIR, "kegg_pathway.tsv"),
        kegg_module  = os.path.join(DICT_DIR, "kegg_module.tsv"),
        kegg_brite   = os.path.join(DICT_DIR, "kegg_brite.tsv"),
        enzyme_dat   = os.path.join(DICT_DIR, "enzyme.dat"),
        pfam_clans   = os.path.join(DICT_DIR, "Pfam-A.clans.tsv"),
        versions     = os.path.join(DICT_DIR, "dictionary_versions.tsv"),
    log:
        os.path.join(DICT_DIR, "download_annotation_dictionaries.log")
    container:
        config["containers"]["conduitr"]
    params:
        dict_dir       = DICT_DIR,
        go_obo_url     = GO_OBO_URL,
        kegg_rest_url  = KEGG_REST_URL,
        enzyme_dat_url = ENZYME_DAT_URL,
        pfam_clans_url = PFAM_CLANS_URL,
    shell:
        r"""
        mkdir -p {params.dict_dir}
        : > {log}

        curl --fail -L -o {output.go_obo}       {params.go_obo_url}                  &>> {log}
        curl --fail -L -o {output.kegg_ko}      {params.kegg_rest_url}/list/ko       &>> {log}
        curl --fail -L -o {output.kegg_pathway} {params.kegg_rest_url}/list/pathway  &>> {log}
        curl --fail -L -o {output.kegg_module}  {params.kegg_rest_url}/list/module   &>> {log}
        curl --fail -L -o {output.kegg_brite}   {params.kegg_rest_url}/list/brite    &>> {log}
        curl --fail -L -o {output.enzyme_dat}   {params.enzyme_dat_url}              &>> {log}

        curl --fail -L -o {output.pfam_clans}.gz {params.pfam_clans_url}             &>> {log}
        gunzip -f {output.pfam_clans}.gz 2>> {log}

        # Record the resolved source + release/version of each dictionary so the
        # conduit object can document where every description came from.
        GO_VER=$(grep -m1 '^data-version:' {output.go_obo} | sed 's/data-version: *//')
        # KEGG's /info/kegg no longer emits a "Release" line, so grep can match
        # nothing; tolerate that (|| true) and fall back to "snapshot" below
        # rather than aborting the rule under set -euo pipefail.
        KEGG_VER=$(curl --fail -sL {params.kegg_rest_url}/info/kegg 2>> {log} \
                   | grep -i 'release' | head -1 | sed 's/^[[:space:]]*//' || true)
        {{
          printf 'dictionary\tsource\tversion\n'
          printf 'go\t%s\t%s\n'     "{params.go_obo_url}"               "$GO_VER"
          printf 'kegg\t%s\t%s\n'   "{params.kegg_rest_url}/list"       "${{KEGG_VER:-snapshot}}"
          printf 'enzyme\t%s\t%s\n' "{params.enzyme_dat_url}"           "snapshot"
          printf 'pfam\t%s\t%s\n'   "{params.pfam_clans_url}"           "Pfam37.0"
        }} > {output.versions}
        """

# Consolidating annotations
rule consolidate_annotations:
  input:
    # Term-name dictionaries (authoritative, for description backfill)
    go_obo = os.path.join(DICT_DIR, "go-basic.obo"),
    kegg_ko = os.path.join(DICT_DIR, "kegg_ko.tsv"),
    kegg_pathway = os.path.join(DICT_DIR, "kegg_pathway.tsv"),
    kegg_module = os.path.join(DICT_DIR, "kegg_module.tsv"),
    kegg_brite = os.path.join(DICT_DIR, "kegg_brite.tsv"),
    enzyme_dat = os.path.join(DICT_DIR, "enzyme.dat"),
    pfam_clans = os.path.join(DICT_DIR, "Pfam-A.clans.tsv"),
    dictionary_versions = os.path.join(DICT_DIR, "dictionary_versions.tsv"),
    # Annotations
    go_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/go_info.txt"),
    kegg_pathway_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/kegg_pathway_info.txt"),
    kegg_map_pathway_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/kegg_map_pathway_info.txt"),
    kegg_orthology_info =  os.path.join(RUN_DIR,"database_resources/detected_protein_resources/kegg_orthology_info.txt"),
    pfam_info =  os.path.join(RUN_DIR,"database_resources/detected_protein_resources/pfam_info.txt"),
    cazy_class_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/cazy_class_info.txt"),
    cazy_family_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/cazy_family_info.txt"),
    eggnog_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/eggnog_info.txt"),
    eggnog_code_info = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/eggnog_code_info.txt"),
    emapper_annotations = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/emapper_annotations.txt"),
    # QF
    qf = os.path.join(RUN_DIR,"output_files/qf.rds")
  output:
    conduit_annotations = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/conduit_annotations.txt")
  log: os.path.join(RUN_DIR,"logs/annotation/consolidate_annotations.log")
  container: config["containers"]["conduitr"]
  script:
    "scripts/consolidate_annotations.R"
