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
taxid_family_map_fp      <- snakemake@input[["taxid_family_map"]]
ncbi_taxonomy_id_fp      <- snakemake@output[["ncbi_taxonomy_id"]]
fdr_results_fp           <- snakemake@output[["fdr_results"]]

conduitR::log_with_timestamp("Input parquet:      %s", first_pass_diann_parquet)
conduitR::log_with_timestamp("taxid→family map:   %s", taxid_family_map_fp)
conduitR::log_with_timestamp("Output (detected):  %s", ncbi_taxonomy_id_fp)
conduitR::log_with_timestamp("Output (FDR table): %s", fdr_results_fp)

# =============================================================================
# Read first-pass DIA-NN results
# =============================================================================
conduitR::log_with_timestamp("Reading first-pass DIA-NN parquet")
precursors <- arrow::read_parquet(first_pass_diann_parquet)
conduitR::log_with_timestamp("Parquet rows: %d", nrow(precursors))

if (nrow(precursors) == 0) {
  stop(sprintf(
    "DIA-NN parquet at %s contains 0 rows — the upstream DIA-NN search produced no peptides. Inspect the corresponding DIA-NN log under logs/ before re-running.",
    first_pass_diann_parquet
  ), call. = FALSE)
}

# =============================================================================
# Extract family taxid and build PSM table for FDR
# =============================================================================
# DIA-NN reports Protein.Ids as "umgap|{id}|{lca_taxid}". We extract the
# lca_taxid and join with taxid_to_family_genus.tsv (col 1 = lca_taxid,
# col 3 = family_taxid) to recover the family for every PSM.
conduitR::log_with_timestamp("Loading taxid → family mapping")
taxid_map <- readr::read_tsv(
  taxid_family_map_fp,
  col_names = c("lca_taxid", "rank", "family_taxid", "genus_taxid"),
  col_types = "cccc",
  show_col_types = FALSE
)
conduitR::log_with_timestamp("taxid map rows: %d", nrow(taxid_map))

conduitR::log_with_timestamp("Extracting lca_taxid from Protein.Ids and mapping to family")

psms <- precursors |>
  dplyr::select(PEP, Stripped.Sequence, Decoy, Protein.Ids) |>
  dplyr::mutate(
    lca_taxid = stringr::str_extract(Protein.Ids, "(?<=\\|)[^|]+$"),
    decoy     = as.logical(Decoy)
  ) |>
  dplyr::left_join(
    dplyr::select(taxid_map, lca_taxid, family_taxid),
    by = "lca_taxid"
  ) |>
  dplyr::filter(!is.na(family_taxid), !is.na(PEP))

conduitR::log_with_timestamp("PSMs with family taxid: %d (targets: %d, decoys: %d)",
  nrow(psms),
  sum(!psms$decoy),
  sum(psms$decoy)
)

if (nrow(psms) == 0) {
  conduitR::log_with_timestamp("WARNING: No PSMs mapped to a family taxid in first-pass results")
}

# =============================================================================
# Apply taxonomic FDR via target-decoy competition at family level
# =============================================================================
conduitR::log_with_timestamp("Running calc_taxon_fdr at family level (FDR threshold = 0.01)")

fdr_result <- conduitR::calc_taxon_fdr(
  pep           = psms$PEP,
  taxon         = psms$family_taxid,
  decoy         = psms$decoy,
  peptide       = psms$Stripped.Sequence,
  fdr_threshold = 0.01
)

conduitR::log_with_timestamp(
  "FDR result: %d target families, %d decoy families, %d detected at FDR<=0.01",
  fdr_result$n_targets, fdr_result$n_decoys, nrow(fdr_result$detected)
)

# =============================================================================
# Format detected families for output
# =============================================================================
detected_families <- fdr_result$detected |>
  dplyr::select(ncbi_taxonomy_id = taxon) |>
  dplyr::mutate(detected_taxonomy = paste0("family_", ncbi_taxonomy_id))

conduitR::log_with_timestamp("Detected %d families", nrow(detected_families))

# =============================================================================
# Write output
# =============================================================================
readr::write_tsv(detected_families, ncbi_taxonomy_id_fp)
conduitR::log_with_timestamp("Written: %s", ncbi_taxonomy_id_fp)

readr::write_tsv(fdr_result$results, fdr_results_fp)
conduitR::log_with_timestamp("Written FDR results table: %s", fdr_results_fp)

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
