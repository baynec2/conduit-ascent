import os

RUN_DIR = config["run_dir"]
EGGNOG_DB_DIR = config["eggnogmapper_db_dir"]

# Required eggNOG-mapper database files
REQUIRED_EGGNOG_FILES = (
    "eggnog.db",
    "eggnog_proteins.dmnd",
)

# Where the emapper database is fetched from.
#
# eggnog6.embl.de 301-redirects to eggnogdb.org, which 404s for both
# emapperdb-5.0.2 files — so the old host does not serve this data at all any
# more, at any path. eggnog5.embl.de serves it directly, with range requests,
# which is what makes the resumed download below work.
#
# Config-overridable so the next host move is a config change rather than a
# code one, matching how external_annotations.smk treats its own URLs.
EGGNOGMAPPER_DB_URL = config.get(
    "eggnogmapper_db_url",
    "http://eggnog5.embl.de/download/emapperdb-5.0.2",
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
        config["containers"]["eggnogmapper"]
    params:
        db_dir = lambda w, output: os.path.dirname(output.db_files[0]),
        db_url = EGGNOGMAPPER_DB_URL,
        # Driven off the same tuple as `output`, so the loop below cannot
        # fetch a different set of files than the rule promises.
        db_files = " ".join(REQUIRED_EGGNOG_FILES)
    shell:
        # This rule used to report success on total failure, and the reason
        # is worth stating precisely, because the obvious diagnosis is wrong.
        #
        # `set -e` was never missing: snakemake already prefixes every shell
        # directive with `set -euo pipefail` (shell.py). The problem is that
        # `set -e` is *ignored* for a failing command in an AND-OR list, so
        # the old `wget -q ... && gunzip ...` swallowed wget's failure
        # regardless:
        #
        #     $ bash -euo pipefail -c 'false && echo b; echo REACHED; exit 0'
        #     REACHED
        #
        # Execution fell through to the final echo and the shell exited 0.
        # The only thing that caught it was snakemake noticing the declared
        # outputs were missing — and `-q` had meanwhile swallowed wget's
        # error, so the log read "Starting..." / "finished!" beside a 0-byte
        # .gz with nothing naming the cause.
        #
        # The fix is therefore the restructuring below — one command per
        # statement, no `&&` chain — not the `set -euo pipefail` line, which
        # duplicates snakemake's and is kept only to make the intent explicit
        # if that prefix ever changes.
        #
        # -nv rather than -q: quiet enough for an 11 GB download, loud enough
        # to record an HTTP error. -c to resume, because these two files are
        # ~6.3 GB and ~4.9 GB and a dropped connection should not restart
        # from zero.
        """
        set -euo pipefail
        mkdir -p {params.db_dir}
        mkdir -p $(dirname {log})
        exec >> {log} 2>&1

        echo "Starting eggNOG-mapper DB download from {params.db_url}..."
        for f in {params.db_files}; do
            echo "Fetching $f.gz"
            wget -nv -c -P {params.db_dir} {params.db_url}/$f.gz
            gunzip -f {params.db_dir}/$f.gz
        done
        echo "eggNOG-mapper DB download finished!"
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
    benchmark:
        os.path.join(RUN_DIR, "benchmarks/annotation/eggnogmapper/run_eggnogmapper.tsv")
    log:
        os.path.join(RUN_DIR, "logs/annotation/eggnogmapper/run_eggnogmapper.log")
    container:
        config["containers"]["eggnogmapper"]
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
        config["containers"]["conduitr"]
    script:
        "scripts/parse_eggnogmapper_annotations.R"
