# =============================================================================
# Setup and Logging
# =============================================================================
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz, append = TRUE)
sink(zz, type = "message")

start_time <- Sys.time()

conduitR::log_with_timestamp("Starting infer_family_presence.R")

# Pure-function helpers (normalize_param, extract_family_psms) live in a
# sibling file so they're reachable from testthat without snakemake@ globals.
snakemake@source("presence_lib.R")

# =============================================================================
# Inputs / Outputs / Config
# =============================================================================
first_pass_diann_parquet <- snakemake@input[["first_pass_diann"]]
taxid_family_map_fp      <- snakemake@input[["taxid_family_map"]]
ncbi_taxonomy_id_fp      <- snakemake@output[["ncbi_taxonomy_id"]]
fdr_results_fp           <- snakemake@output[["fdr_results"]]

# Presence-call params (see config/snakemake.yaml).
# `method` selects the presence rule: "qvalue" (global ranked-list q-value gate)
# or "enrichment" (per-taxon per-peptide score >= margin x pooled decoy noise
# rate — robust at taxon scale; the recommended default).
method                   <- snakemake@params[["method"]]
margin                   <- snakemake@params[["margin"]]
qvalue_threshold         <- snakemake@params[["qvalue_threshold"]]
min_peptides             <- snakemake@params[["min_peptides"]]
min_confident_peptides   <- snakemake@params[["min_confident_peptides"]]
# Optional abundance-based coverage filter (disabled by default).
score_fraction_threshold <- snakemake@params[["score_fraction_threshold"]]
max_taxa                 <- snakemake@params[["max_taxa"]]

if (is.null(method) || length(method) == 0 || is.na(method) || method == "") {
  method <- "enrichment"
}
method <- match.arg(as.character(method), c("qvalue", "enrichment", "count"))
margin                   <- normalize_param(margin)
qvalue_threshold         <- normalize_param(qvalue_threshold)
min_peptides             <- normalize_param(min_peptides)
min_confident_peptides   <- normalize_param(min_confident_peptides)
score_fraction_threshold <- normalize_param(score_fraction_threshold)
max_taxa                 <- normalize_param(max_taxa)

if (is.na(margin))                 margin                 <- 2
if (is.na(qvalue_threshold))       qvalue_threshold       <- 0.05
if (is.na(min_peptides))           min_peptides           <- 2
if (is.na(min_confident_peptides)) min_confident_peptides <- 10

# For method = "count", the picked target-decoy competition is not the gate, but
# we still run it (as "qvalue") to populate the audit table's score/fdr/qvalue
# columns. The count gate is applied in apply_picked_presence_filter.
audit_method <- if (method == "count") "qvalue" else method

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
  "Presence rule: method=%s, margin=%s, qvalue_threshold=%s, min_peptides=%s, min_confident_peptides=%s; coverage filter: score_fraction_threshold=%s, max_taxa=%s",
  method, format(margin), format(qvalue_threshold), format(min_peptides),
  format(min_confident_peptides), format(score_fraction_threshold), format(max_taxa)
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

psms <- extract_family_psms(precursors, taxid_map)

conduitR::log_with_timestamp("PSMs with family taxid: %d (targets: %d, decoys: %d)",
  nrow(psms),
  sum(!psms$decoy),
  sum(psms$decoy)
)

if (nrow(psms) == 0) {
  conduitR::log_with_timestamp("WARNING: No PSMs mapped to a family taxid in first-pass results")
}

# =============================================================================
# Picked target-decoy FDR at family level
# =============================================================================
# conduitR::call_taxon_presence aggregates per-(family, decoy) scores and applies the
# picked target-decoy competition in one call: each family keeps only the
# higher-scoring of its {target, decoy} pair, then all representatives compete
# in one ranked list. This stops a high-abundance family's reversed decoy from
# outranking a true low-abundance family's target. The min_peptides gate is on
# n_unique_peptides_all.
conduitR::log_with_timestamp(
  "Running picked FDR at family level (method=%s, margin=%s, qvalue_threshold=%s, min_peptides=%s)",
  method, format(margin), format(qvalue_threshold), format(min_peptides)
)

fdr_result <- conduitR::call_taxon_presence(
  pep              = psms$PEP,
  taxon            = psms$family_taxid,
  decoy            = psms$decoy,
  peptide          = psms$Stripped.Sequence,
  qvalue_threshold = qvalue_threshold,
  min_peptides     = min_peptides,
  method           = audit_method,
  margin           = margin
)

conduitR::log_with_timestamp(
  "Picked FDR: %d target reps, %d decoy reps, first decoy winner at rank %s of %d",
  fdr_result$n_targets, fdr_result$n_decoys,
  format(fdr_result$first_decoy_rank), fdr_result$n_targets + fdr_result$n_decoys
)
if (fdr_result$n_missing_pair > 0) {
  conduitR::log_with_timestamp(
    "WARNING: %d families present on only one side (target-only or decoy-only); pairing incomplete",
    fdr_result$n_missing_pair
  )
}

# Confident-peptide count per (target) family: distinct peptide sequences each
# with at least one precursor PSM at Q.Value <= 0.01. Informational only —
# carried into the audit table alongside (not as) the presence gate.
PSM_QVALUE_CONFIDENT <- 0.01
q01_peptide_counts <- psms |>
  dplyr::filter(!decoy, !is.na(Q.Value), Q.Value <= PSM_QVALUE_CONFIDENT) |>
  dplyr::distinct(family_taxid, Stripped.Sequence) |>
  dplyr::count(family_taxid, name = "n_unique_peptides_q01") |>
  dplyr::rename(taxon = family_taxid)

# =============================================================================
# Audit table + optional score-coverage filter
# =============================================================================
picked_out <- apply_picked_presence_filter(
  picked_results           = fdr_result$results,
  q01_counts               = q01_peptide_counts,
  min_peptides             = min_peptides,
  score_fraction_threshold = score_fraction_threshold,
  max_taxa                 = max_taxa,
  method                   = method,
  min_confident_peptides   = min_confident_peptides
)

conduitR::log_with_timestamp(
  "first-pass filter: %d families passed picked FDR; %d carried forward (coverage mode=%s, value=%s)",
  picked_out$n_passers, picked_out$n_carried,
  picked_out$filter_mode, format(picked_out$filter_value)
)

augmented_results <- picked_out$augmented

# =============================================================================
# Format detected families for output
# =============================================================================
detected_families <- augmented_results |>
  dplyr::filter(carried_forward) |>
  dplyr::transmute(
    ncbi_taxonomy_id  = taxon,
    detected_taxonomy = paste0("family_", taxon)
  )

conduitR::log_with_timestamp("Carried forward %d families", nrow(detected_families))

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
  # Preserve the historical column order; picked_winner + pass are the new
  # additions (the FDR computation now writes picked fdr / qvalue / pass).
  dplyr::select(
    taxon, taxon_name, taxon_name_rank, score, n_unique_peptides_all, decoy,
    picked_winner, fdr, qvalue, pass, n_unique_peptides_q01,
    score_fraction, cumulative_score_fraction, carried_forward, filter_reason
  )

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
