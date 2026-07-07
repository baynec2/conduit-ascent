################################################################################
# Shared genome-resource cache helpers
################################################################################
# Centralizes the per-(catalog) and per-(catalog,filter,max_genomes) cache-key
# derivation used by genomes.smk, hapid.smk, and genome_peptidotyping.smk so that
# expensive genome-derived artifacts (Prodigal/FGS FAA, HMMER tblout, bakta
# annotations, LCA peptide DBs, hapid speclib) can be reused across runs that
# share the same MGnify reference set.
#
# Sharing only activates when genome_download_source == "mgnify" — MGnify
# accessions are globally unique. For local genome runs the helpers fall back to
# per-run paths so cross-experiment name collisions can't silently mix data.
#
# Include from each consumer module via:
#   include: os.path.join(workflow.basedir, "modules/search_space/_shared/genome_cache.smk")
################################################################################

import os


def _genome_cache_enabled():
    return config.get("genome_download_source") == "mgnify"


def _catalog_slug():
    return config.get("mgnify_catalog", "").replace("/", "_")


def _genome_set_slug():
    """(catalog, taxonomy_filter, max_genomes) → string. Filter and max_genomes
    must be in the slug because MGnify's species_representatives.txt is filtered
    per-run by these knobs; two runs with different filters genuinely have
    different genome sets and must NOT alias to the same cache directory."""
    flt = config.get("mgnify_taxonomy_filter") or "none"
    flt = str(flt).replace("/", "_").replace(" ", "_")
    mx = config.get("mgnify_max_genomes", 0)
    return f"{_catalog_slug()}__filter-{flt}__max{mx}"


def per_genome_cache_root(tool):
    """Root for per-genome artifacts (e.g. one subdir per genome ID).
    Shared across runs when MGnify; per-run otherwise."""
    if _genome_cache_enabled():
        return os.path.join(config["genome_resource_dir"], _catalog_slug(), tool)
    return os.path.join(config["run_dir"], "database_resources", tool)


def per_genome_set_cache_root(module):
    """Root for per-genome-set aggregate artifacts.
    Shared across runs that match the genome-set slug; per-run otherwise."""
    if _genome_cache_enabled():
        return os.path.join(config["genome_set_resource_dir"],
                            _genome_set_slug(), module)
    return os.path.join(config["run_dir"], "database_resources", module)
