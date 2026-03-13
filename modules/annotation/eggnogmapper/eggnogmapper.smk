import os

RUN_DIR = config["run_dir"]
EGGNOG_DB_DIR = config["eggnogmapper_db_dir"]

# Required eggNOG-mapper database files
REQUIRED_EGGNOG_FILES = (
    "eggnog.db",
    "eggnog_proteins.dmnd",
)

################################################################################
# Download eggNOG-mapper database (once)
################################################################################
rule download_eggnogmapper_db:
    output:
        db_files = expand(os.path.join(EGGNOG_DB_DIR, "{file}"), file=REQUIRED_EGGNOG_FILES)
    log:
        os.path.join(EGGNOG_DB_DIR, "logs/download_eggnogmapper_db.log")
    container:
        "docker://baynec2/eggnogmapper:2.1.12"
    params:
        db_dir = lambda w, output: os.path.dirname(output.db_files[0])
    shell:
        """
        mkdir -p {params.db_dir}
        mkdir -p $(dirname {log})
        echo "Starting eggNOG-mapper DB download..." > {log}

        wget -q -O {params.db_dir}/eggnog.db.gz \
            http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.db.gz >> {log} 2>&1 \
            && gunzip {params.db_dir}/eggnog.db.gz >> {log} 2>&1

        wget -q -O {params.db_dir}/eggnog_proteins.dmnd.gz \
            http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz >> {log} 2>&1 \
            && gunzip {params.db_dir}/eggnog_proteins.dmnd.gz >> {log} 2>&1

        echo "eggNOG-mapper DB download finished!" >> {log}
        """

################################################################################
# Run eggNOG-mapper on detected proteins
################################################################################
rule run_eggnogmapper:
    input:
        fasta = os.path.join(RUN_DIR, "database_resources/detected_protein_resources/detected_protein.fasta"),
        db_files = expand(os.path.join(EGGNOG_DB_DIR, "{file}"), file=REQUIRED_EGGNOG_FILES)
    output:
        annotations = os.path.join(RUN_DIR, "database_resources/detected_protein_resources/emapper.emapper.annotations")
    log:
        os.path.join(RUN_DIR, "logs/annotation/eggnogmapper/run_eggnogmapper.log")
    container:
        "docker://baynec2/eggnogmapper:2.1.12"
    params:
        db_dir = lambda w, input: os.path.dirname(input.db_files[0]),
        output_dir = lambda w, output: os.path.dirname(output.annotations),
        output_prefix = "emapper"
    threads: workflow.cores
    shell:
        """
        mkdir -p $(dirname {log})
        echo "Running eggNOG-mapper..." > {log}

        emapper.py \
            -i {input.fasta} \
            --itype proteins \
            --output {params.output_prefix} \
            --output_dir {params.output_dir} \
            --data_dir {params.db_dir} \
            --cpu {threads} \
            --override >> {log} 2>&1

        echo "eggNOG-mapper finished!" >> {log}
        """

################################################################################
# Parse eggNOG-mapper output into standard long-format annotation table
################################################################################
rule parse_eggnogmapper_annotations:
    input:
        annotations = os.path.join(RUN_DIR, "database_resources/detected_protein_resources/emapper.emapper.annotations")
    output:
        emapper_annotations = os.path.join(RUN_DIR, "database_resources/detected_protein_resources/emapper_annotations.txt")
    log:
        os.path.join(RUN_DIR, "logs/annotation/eggnogmapper/parse_eggnogmapper_annotations.log")
    container:
        "docker://baynec2/conduitr:alpha"
    script:
        "scripts/parse_eggnogmapper_annotations.R"
