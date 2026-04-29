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

# Score-coverage filter params (see config/snakemake.yaml).
score_fraction_threshold <- snakemake@params[["score_fraction_threshold"]]
max_taxa                 <- snakemake@params[["max_taxa"]]

normalize_param <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  if (is.character(x) && (x %in% c("", "NA", "null", "None"))) return(NA_real_)
  suppressWarnings(as.numeric(x))
}
score_fraction_threshold <- normalize_param(score_fraction_threshold)
max_taxa                 <- normalize_param(max_taxa)

if (!is.na(score_fraction_threshold) && !is.na(max_taxa)) {
  stop(sprintf(
    paste0(
      "Both peptidotyping_second_pass_score_fraction_threshold (%s) and ",
      "peptidotyping_second_pass_max_taxa (%s) are set in the config. ",
      "Set exactly one per pass (or both null to disable)."
    ),
    format(score_fraction_threshold), format(max_taxa)
  ))
}

conduitR::log_with_timestamp("Input parquet:      %s", second_pass_diann_parquet)
conduitR::log_with_timestamp("Output (detected):  %s", detected_species_strains_fp)
conduitR::log_with_timestamp("Output (FDR table): %s", fdr_results_fp)
conduitR::log_with_timestamp(
  "Score-coverage filter: score_fraction_threshold=%s, max_taxa=%s",
  format(score_fraction_threshold), format(max_taxa)
)

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
# Apply score-coverage filter to FDR-passing species/strains
# =============================================================================
# Sort FDR-passers by score desc, compute per-taxon score fractions, and decide
# which taxa are carried forward to ncbi_taxa_ids.txt per the configured filter.
candidates <- fdr_result$results |>
  dplyr::filter(!decoy, qvalue <= 0.01) |>
  dplyr::arrange(dplyr::desc(score))

n_candidates <- nrow(candidates)

if (n_candidates == 0) {
  candidates <- candidates |>
    dplyr::mutate(
      score_fraction            = numeric(0),
      cumulative_score_fraction = numeric(0),
      carried_forward           = logical(0)
    )
  filter_mode  <- "none (no FDR-passing taxa)"
  filter_value <- NA_real_
} else {
  total_score <- sum(candidates$score)
  candidates$score_fraction            <- candidates$score / total_score
  candidates$cumulative_score_fraction <- cumsum(candidates$score_fraction)

  if (!is.na(max_taxa)) {
    cutoff       <- min(as.integer(max_taxa), n_candidates)
    filter_mode  <- "max_taxa"
    filter_value <- as.numeric(max_taxa)
  } else if (!is.na(score_fraction_threshold)) {
    reach <- which(candidates$cumulative_score_fraction >= score_fraction_threshold)
    cutoff       <- if (length(reach) == 0) n_candidates else reach[1]
    filter_mode  <- "score_fraction_threshold"
    filter_value <- score_fraction_threshold
  } else {
    cutoff       <- n_candidates
    filter_mode  <- "none (filter disabled)"
    filter_value <- NA_real_
  }

  candidates$carried_forward <- seq_len(n_candidates) <= cutoff
}

n_carried <- sum(candidates$carried_forward)
conduitR::log_with_timestamp(
  "second-pass filter: kept %d/%d FDR-passing species/strains (mode=%s, value=%s)",
  n_carried, n_candidates, filter_mode, format(filter_value)
)

# =============================================================================
# Format detected species/strains for output
# =============================================================================
detected_species_strains <- candidates |>
  dplyr::filter(carried_forward) |>
  dplyr::transmute(
    ncbi_taxonomy_id  = taxon,
    detected_taxonomy = paste0("species_strain_", taxon)
  )

conduitR::log_with_timestamp("Carried forward %d species/strains", nrow(detected_species_strains))

# =============================================================================
# Build augmented FDR results table (audit trail)
# =============================================================================
non_candidates <- fdr_result$results |>
  dplyr::filter(decoy | qvalue > 0.01) |>
  dplyr::mutate(
    score_fraction            = NA_real_,
    cumulative_score_fraction = NA_real_,
    carried_forward           = FALSE
  )

augmented_results <- dplyr::bind_rows(candidates, non_candidates)

# =============================================================================
# Write output
# =============================================================================
readr::write_tsv(detected_species_strains, detected_species_strains_fp)
conduitR::log_with_timestamp("Written: %s", detected_species_strains_fp)

readr::write_tsv(augmented_results, fdr_results_fp)
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
