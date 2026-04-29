# =============================================================================
# Setup and Logging
# =============================================================================
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz, append = TRUE)
sink(zz, type = "message")

start_time <- Sys.time()
conduitR::log_with_timestamp("Starting infer_species_presence.R")

# =============================================================================
# Inputs / Outputs / Config
# =============================================================================
hapid_diann_parquet  <- snakemake@input[["hapid_diann"]]
ncbi_taxa_ids_fp     <- snakemake@output[["ncbi_taxa_ids"]]
presence_min_peptides <- snakemake@config[["presence_min_peptides"]]

conduitR::log_with_timestamp("Input parquet: %s", hapid_diann_parquet)
conduitR::log_with_timestamp("Output: %s", ncbi_taxa_ids_fp)
conduitR::log_with_timestamp("presence_min_peptides threshold: %d", presence_min_peptides)

# =============================================================================
# Read DIA-NN results and infer species/strain presence
# =============================================================================
conduitR::log_with_timestamp("Reading HAPiID first-pass DIA-NN parquet")
precursors <- arrow::read_parquet(hapid_diann_parquet)

# UMGAP first-pass DB headers have format `umgap|{id}|{lca_il}` — DIA-NN
# surfaces this as Protein.Group, with the trailing `|<lca_il>` carrying the
# species/strain NCBI taxon ID. Protein.Names is the human-readable label
# (`species_<name>`) and does NOT contain an OX= token, so don't read from it.
detected <- precursors |>
  dplyr::filter(Proteotypic == 1) |>
  dplyr::mutate(
    ncbi_taxonomy_id = sub(".*\\|", "", Protein.Group)
  ) |>
  dplyr::group_by(ncbi_taxonomy_id, Stripped.Sequence) |>
  dplyr::summarise(.groups = "drop") |>
  dplyr::group_by(ncbi_taxonomy_id) |>
  dplyr::summarise(n_peptides = dplyr::n(), .groups = "drop") |>
  dplyr::filter(n_peptides >= presence_min_peptides) |>
  dplyr::select(ncbi_taxonomy_id) |>
  dplyr::distinct()

conduitR::log_with_timestamp("Detected %d species/strains above threshold", nrow(detected))

# Write ncbi_taxa_ids.txt — consumed by ncbi_taxonomy_id workflow
# Canonical single-column format with header `ncbi_taxonomy_id`.
readr::write_delim(detected, ncbi_taxa_ids_fp)

# =============================================================================
# Cleanup and Logging
# =============================================================================
end_time <- Sys.time()
elapsed_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
conduitR::log_with_timestamp(
  "Completed infer_species_presence.R. Time taken: %.2f minutes", elapsed_minutes
)

sink(type = "message")
sink()
close(zz)
