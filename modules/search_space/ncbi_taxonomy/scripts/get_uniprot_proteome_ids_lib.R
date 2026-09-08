# Pure-function helper for get_uniprot_proteome_ids.R. Kept in a separate file
# so the TSV-ingest + append logic is reachable from testthat without
# snakemake@ globals.

# Parse organism IDs from an already-read TSV tibble. Tolerates either a
# single-column format (just the numeric ID under any header, e.g.
# `ncbi_taxa_id`) or a multi-column TSV where the first column carries the IDs
# and additional columns carry user metadata (e.g. `organism_id\tsource`).
# The previous file-level implementation used `readr::read_lines + as.integer`,
# which silently coerced lines like `820\tBacteroides_uniformis_ATCC_8492` to
# NA and then queried UniProt with NAs (HTTP 400).
#
# `append_id` is the value from snakemake@config$append_additional_ncbi_taxa_id
# — either `FALSE` (default; nothing to append) or a numeric/integer ID.
#
# Returns a deduplicated integer vector of organism IDs.
parse_organism_ids <- function(raw_df, append_id = FALSE) {
  organism_ids <- raw_df |>
    dplyr::pull(1) |>
    (\(x) suppressWarnings(as.integer(x)))() |>
    (\(x) x[!is.na(x)])() |>
    unique()

  if (!isFALSE(append_id) && !(append_id %in% organism_ids)) {
    organism_ids <- c(organism_ids, as.integer(append_id))
  }

  organism_ids
}
