# =============================================================================
# Setup and Logging
# =============================================================================
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz, append = TRUE)
sink(zz, type = "message")

start_time <- Sys.time()
conduitR::log_with_timestamp("Starting infer_family_presence.R")

# =============================================================================
# Inputs / Outputs / Config
# =============================================================================
first_pass_diann_parquet <- snakemake@input[["first_pass_diann"]]
rank_mapping_fp          <- snakemake@input[["rank_mapping"]]
ncbi_taxonomy_id_fp      <- snakemake@output[["ncbi_taxonomy_id"]]
presence_min_peptides    <- snakemake@config[["presence_min_peptides"]]

conduitR::log_with_timestamp("Input parquet:  %s", first_pass_diann_parquet)
conduitR::log_with_timestamp("Rank mapping:   %s", rank_mapping_fp)
conduitR::log_with_timestamp("Output:         %s", ncbi_taxonomy_id_fp)
conduitR::log_with_timestamp("presence_min_peptides: %d", presence_min_peptides)

# =============================================================================
# Read effective detection rank mapping
# =============================================================================
# Columns: family_taxid, effective_rank, n_peptides, representative_taxid
conduitR::log_with_timestamp("Reading effective detection rank mapping")
rank_mapping <- readr::read_tsv(rank_mapping_fp, col_types = readr::cols(.default = "c"))
conduitR::log_with_timestamp(
  "Rank mapping: %d families (%d family, %d genus, %d species_strain)",
  nrow(rank_mapping),
  sum(rank_mapping$effective_rank == "family"),
  sum(rank_mapping$effective_rank == "genus"),
  sum(rank_mapping$effective_rank == "species_strain")
)

# =============================================================================
# Read first-pass DIA-NN results
# =============================================================================
conduitR::log_with_timestamp("Reading first-pass DIA-NN parquet: %s", first_pass_diann_parquet)
precursors <- arrow::read_parquet(first_pass_diann_parquet)

# =============================================================================
# Extract FAM= family taxid from every FASTA header hit
# =============================================================================
# All entries in effective_first_pass_database.fasta have a FAM=<taxid> tag
# appended to their header. DIA-NN reports these headers in Protein.Names.
# We extract FAM= to always recover the family taxid, regardless of whether
# the detection was at family, genus, or species/strain level.
conduitR::log_with_timestamp("Extracting FAM= family taxid from Protein.Names")

pep <- precursors |>
  dplyr::filter(Proteotypic == 1) |>
  dplyr::mutate(
    family_taxid = stringr::str_extract(Protein.Names, "(?<=FAM=)[0-9]+"),
    detected_taxonomy = gsub(".*_", "", Protein.Names)
  ) |>
  dplyr::filter(!is.na(family_taxid))

if (nrow(pep) == 0) {
  conduitR::log_with_timestamp("WARNING: No proteotypic hits with FAM= tag found in first-pass results")
}

# =============================================================================
# Count unique peptides per family across all runs
# =============================================================================
conduitR::log_with_timestamp("Counting unique peptides per family")

detected_families <- pep |>
  dplyr::group_by(family_taxid, Stripped.Sequence) |>
  dplyr::summarise(sum_intensity = sum(Precursor.Normalised), .groups = "drop") |>
  dplyr::group_by(family_taxid) |>
  dplyr::summarise(n_peptides = dplyr::n(), .groups = "drop") |>
  dplyr::filter(n_peptides >= presence_min_peptides) |>
  dplyr::rename(ncbi_taxonomy_id = family_taxid) |>
  dplyr::mutate(detected_taxonomy = paste0("family_", ncbi_taxonomy_id)) |>
  dplyr::select(ncbi_taxonomy_id, detected_taxonomy) |>
  dplyr::distinct()

conduitR::log_with_timestamp("Detected %d families above threshold (%d peptides)",
  nrow(detected_families), presence_min_peptides)

# Log breakdown by effective rank for detected families
detected_with_rank <- detected_families |>
  dplyr::left_join(
    rank_mapping |> dplyr::select(family_taxid, effective_rank) |>
      dplyr::rename(ncbi_taxonomy_id = family_taxid),
    by = "ncbi_taxonomy_id"
  )

conduitR::log_with_timestamp(
  "Detection breakdown: %d via family peptides, %d via genus peptides, %d via species/strain peptides",
  sum(detected_with_rank$effective_rank == "family", na.rm = TRUE),
  sum(detected_with_rank$effective_rank == "genus", na.rm = TRUE),
  sum(detected_with_rank$effective_rank == "species_strain", na.rm = TRUE)
)

# =============================================================================
# Write output
# =============================================================================
readr::write_delim(detected_families, ncbi_taxonomy_id_fp)
conduitR::log_with_timestamp("Written: %s", ncbi_taxonomy_id_fp)

# =============================================================================
# Cleanup and Logging
# =============================================================================
end_time <- Sys.time()
elapsed_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
conduitR::log_with_timestamp(
  "Completed infer_family_presence.R. Time taken: %.2f minutes", elapsed_minutes
)

sink(type = "message")
sink()
close(zz)
