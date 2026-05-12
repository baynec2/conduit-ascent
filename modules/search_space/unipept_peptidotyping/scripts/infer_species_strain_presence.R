# =============================================================================
# Setup and Logging
# =============================================================================
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz, append = TRUE)
sink(zz, type = "message")

start_time <- Sys.time()

conduitR::log_with_timestamp("Starting infer_species_strain_presence.R")

# =============================================================================
# Inputs / Outputs
# =============================================================================
second_pass_diann_parquet   <- snakemake@input[["second_pass_diann"]]
detected_species_strains_fp <- snakemake@output[["detected_species_strains"]]
fdr_results_fp              <- snakemake@output[["fdr_results"]]

conduitR::log_with_timestamp("Input parquet:      %s", second_pass_diann_parquet)
conduitR::log_with_timestamp("Output (detected):  %s", detected_species_strains_fp)
conduitR::log_with_timestamp("Output (FDR table): %s", fdr_results_fp)

# =============================================================================
# Read second-pass DIA-NN results
# =============================================================================
conduitR::log_with_timestamp("Reading second-pass DIA-NN parquet")
precursors <- arrow::read_parquet(second_pass_diann_parquet)
conduitR::log_with_timestamp("Parquet rows: %d", nrow(precursors))

# =============================================================================
# Extract species/strain taxid and build PSM table for FDR
# =============================================================================
# DIA-NN reports Protein.Ids as "umgap|{id}|{lca_taxid}". For second-pass
# entries the lca_taxid IS the species/strain taxid, so we extract it directly.
conduitR::log_with_timestamp("Extracting species/strain taxid from Protein.Ids")

psms <- precursors |>
  dplyr::select(PEP, Stripped.Sequence, Decoy, Protein.Ids) |>
  dplyr::mutate(
    species_taxid = stringr::str_extract(Protein.Ids, "(?<=\\|)[^|]+$"),
    decoy         = as.logical(Decoy)
  ) |>
  dplyr::filter(!is.na(species_taxid), !is.na(PEP))

conduitR::log_with_timestamp("PSMs with species/strain taxid: %d (targets: %d, decoys: %d)",
  nrow(psms),
  sum(!psms$decoy),
  sum(psms$decoy)
)

if (nrow(psms) == 0) {
  conduitR::log_with_timestamp("WARNING: No PSMs with species/strain taxid found in second-pass results")
}

# =============================================================================
# Apply taxonomic FDR via target-decoy competition at species/strain level
# =============================================================================
conduitR::log_with_timestamp("Running calc_taxon_fdr at species/strain level (FDR threshold = 0.01)")

fdr_result <- conduitR::calc_taxon_fdr(
  pep           = psms$PEP,
  taxon         = psms$species_taxid,
  decoy         = psms$decoy,
  peptide       = psms$Stripped.Sequence,
  fdr_threshold = 0.01
)

conduitR::log_with_timestamp(
  "FDR result: %d target taxa, %d decoy taxa, %d detected at FDR<=0.01",
  fdr_result$n_targets, fdr_result$n_decoys, nrow(fdr_result$detected)
)

# =============================================================================
# Format detected species/strains for output
# =============================================================================
detected_species_strains <- fdr_result$detected |>
  dplyr::select(ncbi_taxonomy_id = taxon) |>
  dplyr::mutate(detected_taxonomy = paste0("species_strain_", ncbi_taxonomy_id))

conduitR::log_with_timestamp("Detected %d species/strains", nrow(detected_species_strains))

# =============================================================================
# Write output
# =============================================================================
readr::write_tsv(detected_species_strains, detected_species_strains_fp)
conduitR::log_with_timestamp("Written: %s", detected_species_strains_fp)

readr::write_tsv(fdr_result$results, fdr_results_fp)
conduitR::log_with_timestamp("Written FDR results table: %s", fdr_results_fp)

# =============================================================================
# Cleanup and Logging
# =============================================================================
end_time <- Sys.time()
elapsed_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
conduitR::log_with_timestamp(
  "Completed infer_species_strain_presence.R. Time taken: %.2f minutes", elapsed_minutes
)

sink(type = "message")
sink()
close(zz)
