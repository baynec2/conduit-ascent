################################################################################
# MetaPhlAn Database Management
################################################################################
rule download_metaphlan_resources:
    """
    Download and install MetaPhlAn database resources.
    This rule downloads the latest MetaPhlAn database and creates necessary directory structure.
    """
    output:
        database_dir = directory("resources/metaphlan_databases"),
    container: "containers/metaphlan.def"
    log: "resources/metaphlan_databases/logs/download_metaphlan_resources.log"
    #container: "docker://gmtscience/metaphlan4"
    shell:
        """
        metaphlan --install \
            --bowtie2db {output.database_dir} \
            >> {log} 2>&1
        """

