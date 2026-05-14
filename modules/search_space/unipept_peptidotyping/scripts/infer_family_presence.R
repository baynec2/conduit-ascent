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

# Score-coverage filter params (see config/snakemake.yaml).
score_fraction_threshold <- snakemake@params[["score_fraction_threshold"]]
max_taxa                 <- snakemake@params[["max_taxa"]]
# Evidence-quantity floor: minimum distinct peptide sequences supporting a taxon
# call. Applied as a hard gate after FDR and before the score-coverage filter.
min_unique_peptides      <- snakemake@params[["min_unique_peptides"]]

normalize_param <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  if (is.character(x) && (x %in% c("", "NA", "null", "None"))) return(NA_real_)
  suppressWarnings(as.numeric(x))
}
score_fraction_threshold <- normalize_param(score_fraction_threshold)
max_taxa                 <- normalize_param(max_taxa)
min_unique_peptides      <- normalize_param(min_unique_peptides)

if (!is.na(score_fraction_threshold) && !is.na(max_taxa)) {
  stop(sprintf(
    paste0(
      "Both peptidotyping_first_pass_score_fraction_threshold (%s) and ",
      "peptidotyping_first_pass_max_taxa (%s) are set in the config. ",
      "Set exactly one per pass (or both null to disable)."
    ),
    format(score_fraction_threshold), format(max_taxa)
  ))
}

conduitR::log_with_timestamp("Input parquet:      %s", first_pass_diann_parquet)
conduitR::log_with_timestamp("taxid→family map:   %s", taxid_family_map_fp)
conduitR::log_with_timestamp("Output (detected):  %s", ncbi_taxonomy_id_fp)
conduitR::log_with_timestamp("Output (FDR table): %s", fdr_results_fp)
conduitR::log_with_timestamp(
  "Filters: min_unique_peptides=%s, score_fraction_threshold=%s, max_taxa=%s",
  format(min_unique_peptides),
  format(score_fraction_threshold), format(max_taxa)
)

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
  dplyr::select(PEP, Q.Value, Stripped.Sequence, Decoy, Protein.Ids) |>
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

# Confident-peptide count per (target) family: distinct peptide sequences each
# with at least one precursor PSM at Q.Value <= 0.01. See the second-pass
# script for the rationale (gating uses this column, not n_unique_peptides_all).
PSM_QVALUE_CONFIDENT <- 0.01
q01_peptide_counts <- psms |>
  dplyr::filter(!decoy, !is.na(Q.Value), Q.Value <= PSM_QVALUE_CONFIDENT) |>
  dplyr::distinct(family_taxid, Stripped.Sequence) |>
  dplyr::count(family_taxid, name = "n_unique_peptides_q01") |>
  dplyr::rename(taxon = family_taxid)

fdr_result$results <- fdr_result$results |>
  dplyr::left_join(q01_peptide_counts, by = "taxon") |>
  dplyr::mutate(n_unique_peptides_q01 = ifelse(is.na(n_unique_peptides_q01),
                                               0L, as.integer(n_unique_peptides_q01)))

# =============================================================================
# Apply min-unique-peptides + score-coverage filters to FDR-passing families
# =============================================================================
# Filter chain (see infer_species_strain_presence.R for the same structure):
#   1. FDR / decoy        — qvalue <= 0.01, !decoy
#   2. min_unique_peptides — evidence-quantity floor (NA disables)
#   3. score_fraction     — coverage filter on what survives (1)+(2)
fdr_passers <- fdr_result$results |>
  dplyr::filter(!decoy, qvalue <= 0.01) |>
  dplyr::arrange(dplyr::desc(score))

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
  "first-pass filter: kept %d/%d FDR-passing families (min_unique_peptides=%s rejected %d; coverage mode=%s, value=%s)",
  n_carried, nrow(fdr_passers),
  format(min_unique_peptides), nrow(pep_failures),
  filter_mode, format(filter_value)
)

# =============================================================================
# Format detected families for output
# =============================================================================
detected_families <- candidates |>
  dplyr::filter(carried_forward) |>
  dplyr::transmute(
    ncbi_taxonomy_id  = taxon,
    detected_taxonomy = paste0("family_", taxon)
  )

conduitR::log_with_timestamp("Carried forward %d families", nrow(detected_families))

# =============================================================================
# Build augmented FDR results table (audit trail)
# =============================================================================
# Every calc_taxon_fdr row is preserved, tagged with carried_forward + the
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

# Pull a human-readable name for each family taxid. Every peptide carries
# `lca_taxid` (last segment of Protein.Ids) and `{rank}_{name}` (Protein.Names).
# Joining lca_taxid → family_taxid via taxid_map gives us a per-family pool of
# names; we pick one with rank preference family > genus > species > strain so
# families whose effective first-pass rank is genus/species fall back to a
# genus/species/strain name observed for them rather than NA.
rank_priority <- c("family" = 1L, "genus" = 2L, "species" = 3L, "strain" = 4L)

name_lookup <- precursors |>
  dplyr::transmute(
    lca_taxid      = stringr::str_extract(Protein.Ids, "(?<=\\|)[^|]+$"),
    name_rank      = stringr::str_extract(Protein.Names, "^(family|genus|species|strain)(?=_)"),
    name_stripped  = stringr::str_remove(Protein.Names, "^(family_|genus_|species_|strain_)")
  ) |>
  dplyr::filter(
    !is.na(lca_taxid),
    !is.na(name_rank),
    !is.na(name_stripped),
    name_stripped != ""
  ) |>
  dplyr::inner_join(
    dplyr::select(taxid_map, lca_taxid, family_taxid),
    by = "lca_taxid"
  ) |>
  dplyr::mutate(prio = rank_priority[name_rank]) |>
  dplyr::arrange(family_taxid, prio) |>
  dplyr::distinct(family_taxid, .keep_all = TRUE) |>
  dplyr::transmute(
    taxon           = family_taxid,
    taxon_name      = name_stripped,
    taxon_name_rank = name_rank
  )

augmented_results <- augmented_results |>
  dplyr::left_join(name_lookup, by = "taxon") |>
  dplyr::relocate(taxon_name, taxon_name_rank, .after = taxon)

# =============================================================================
# Write output
# =============================================================================
readr::write_tsv(detected_families, ncbi_taxonomy_id_fp)
conduitR::log_with_timestamp("Written: %s", ncbi_taxonomy_id_fp)

readr::write_tsv(augmented_results, fdr_results_fp)
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
