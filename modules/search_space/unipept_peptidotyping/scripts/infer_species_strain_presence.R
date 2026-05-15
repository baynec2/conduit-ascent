# =============================================================================
# Setup and Logging
# =============================================================================
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz, append = TRUE)
sink(zz, type = "message")

start_time <- Sys.time()

conduitR::log_with_timestamp("Starting infer_species_strain_presence.R")

# Pure-function helpers (normalize_param, extract_species_strain_psms) live in
# a sibling file so they're reachable from testthat without snakemake@ globals.
snakemake@source("presence_lib.R")

# =============================================================================
# Inputs / Outputs
# =============================================================================
second_pass_diann_parquet   <- snakemake@input[["second_pass_diann"]]
detected_species_strains_fp <- snakemake@output[["detected_species_strains"]]
fdr_results_fp              <- snakemake@output[["fdr_results"]]

# Score-coverage filter params (see config/snakemake.yaml).
score_fraction_threshold <- snakemake@params[["score_fraction_threshold"]]
max_taxa                 <- snakemake@params[["max_taxa"]]
# Evidence-quantity floor: minimum distinct peptide sequences supporting a taxon
# call. Applied as a hard gate after FDR and before the score-coverage filter.
min_unique_peptides      <- snakemake@params[["min_unique_peptides"]]

score_fraction_threshold <- normalize_param(score_fraction_threshold)
max_taxa                 <- normalize_param(max_taxa)
min_unique_peptides      <- normalize_param(min_unique_peptides)

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
  "Filters: min_unique_peptides=%s, score_fraction_threshold=%s, max_taxa=%s",
  format(min_unique_peptides),
  format(score_fraction_threshold), format(max_taxa)
)

# =============================================================================
# Read second-pass DIA-NN results
# =============================================================================
conduitR::log_with_timestamp("Reading second-pass DIA-NN parquet")
precursors <- arrow::read_parquet(second_pass_diann_parquet)
conduitR::log_with_timestamp("Parquet rows: %d", nrow(precursors))

if (nrow(precursors) == 0) {
  stop(sprintf(
    "DIA-NN parquet at %s contains 0 rows — the upstream DIA-NN search produced no peptides. Inspect the corresponding DIA-NN log under logs/ before re-running.",
    second_pass_diann_parquet
  ), call. = FALSE)
}

# =============================================================================
# Extract species/strain taxid and build PSM table for FDR
# =============================================================================
# DIA-NN reports Protein.Ids as "umgap|{id}|{lca_taxid}". For second-pass
# entries the lca_taxid IS the species/strain taxid, so we extract it directly.
conduitR::log_with_timestamp("Extracting species/strain taxid from Protein.Ids")

psms <- extract_species_strain_psms(precursors)

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

# Confident-peptide count per (target) taxon: number of distinct peptide
# sequences each with at least one precursor PSM at Q.Value <= 0.01. This is
# the count users typically mean when they say "N peptides per taxon"; it sits
# alongside the broader n_unique_peptides_all from calc_taxon_fdr (which
# includes weak-PEP contributions) so the audit trail surfaces both.
PSM_QVALUE_CONFIDENT <- 0.01
q01_peptide_counts <- psms |>
  dplyr::filter(!decoy, !is.na(Q.Value), Q.Value <= PSM_QVALUE_CONFIDENT) |>
  dplyr::distinct(species_taxid, Stripped.Sequence) |>
  dplyr::count(species_taxid, name = "n_unique_peptides_q01") |>
  dplyr::rename(taxon = species_taxid)

fdr_result$results <- fdr_result$results |>
  dplyr::left_join(q01_peptide_counts, by = "taxon") |>
  dplyr::mutate(n_unique_peptides_q01 = ifelse(is.na(n_unique_peptides_q01),
                                               0L, as.integer(n_unique_peptides_q01)))

# =============================================================================
# Apply min-unique-peptides + score-coverage filters to FDR-passing taxa
# =============================================================================
# Filter chain (each rejection reason is recorded in `filter_reason`):
#   1. FDR / decoy        — handled by upstream filter (qvalue <= 0.01, !decoy)
#   2. min_unique_peptides — evidence-quantity floor (NA disables)
#   3. score_fraction     — coverage filter on what survives (1)+(2)
fdr_passers <- fdr_result$results |>
  dplyr::filter(!decoy, qvalue <= 0.01) |>
  dplyr::arrange(dplyr::desc(score))

# Apply min_unique_peptides as a hard floor — gates on the confident-peptide
# count (n_unique_peptides_q01), not n_unique_peptides_all. Reasoning: the
# filter exists to require real evidence, and a peptide that doesn't clear its
# own precursor q-value isn't standalone evidence.
if (!is.na(min_unique_peptides)) {
  candidates    <- dplyr::filter(fdr_passers, n_unique_peptides_q01 >= min_unique_peptides)
  pep_failures  <- dplyr::filter(fdr_passers, n_unique_peptides_q01 <  min_unique_peptides)
} else {
  candidates    <- fdr_passers
  pep_failures  <- fdr_passers[0, , drop = FALSE]
}

n_candidates <- nrow(candidates)

if (n_candidates == 0) {
  candidates <- candidates |>
    dplyr::mutate(
      score_fraction            = numeric(0),
      cumulative_score_fraction = numeric(0),
      carried_forward           = logical(0),
      filter_reason             = character(0)
    )
  filter_mode  <- "none (no FDR+min-peptides-passing taxa)"
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

  in_cov <- seq_len(n_candidates) <= cutoff
  candidates$carried_forward <- in_cov
  candidates$filter_reason   <- ifelse(in_cov, "", filter_mode)
}

# Tag min-peptide-rejected taxa with NA fractions and explicit reason.
if (nrow(pep_failures) > 0) {
  pep_failures <- pep_failures |>
    dplyr::mutate(
      score_fraction            = NA_real_,
      cumulative_score_fraction = NA_real_,
      carried_forward           = FALSE,
      filter_reason             = "min_unique_peptides"
    )
}

n_carried <- sum(candidates$carried_forward)
conduitR::log_with_timestamp(
  "second-pass filter: kept %d/%d FDR-passing species/strains (min_unique_peptides=%s rejected %d; coverage mode=%s, value=%s)",
  n_carried, nrow(fdr_passers),
  format(min_unique_peptides), nrow(pep_failures),
  filter_mode, format(filter_value)
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
# All calc_taxon_fdr rows are preserved, tagged with carried_forward + the
# filter_reason that explains why a row was rejected (decoy / fdr /
# min_unique_peptides / score_fraction_threshold / max_taxa). Empty reason
# means the taxon passed every gate.
non_candidates <- fdr_result$results |>
  dplyr::filter(decoy | qvalue > 0.01) |>
  dplyr::mutate(
    score_fraction            = NA_real_,
    cumulative_score_fraction = NA_real_,
    carried_forward           = FALSE,
    filter_reason             = ifelse(decoy, "decoy", "fdr")
  )

augmented_results <- dplyr::bind_rows(candidates, pep_failures, non_candidates)

# Pull a human-readable name for each taxon from DIA-NN's Protein.Names column
# (sourced from the FASTA "{rank}_{name}" description field). Taxa with no
# parquet match (e.g. detected families whose FASTA entries were genus/species
# fallbacks only) get NA.
name_map <- precursors |>
  dplyr::transmute(
    taxon      = stringr::str_extract(Protein.Ids, "(?<=\\|)[^|]+$"),
    taxon_name = stringr::str_remove(Protein.Names, "^(family_|genus_|species_|strain_)")
  ) |>
  dplyr::filter(!is.na(taxon), !is.na(taxon_name), taxon_name != "") |>
  dplyr::distinct(taxon, .keep_all = TRUE)

augmented_results <- augmented_results |>
  dplyr::left_join(name_map, by = "taxon") |>
  dplyr::relocate(taxon_name, .after = taxon)

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
