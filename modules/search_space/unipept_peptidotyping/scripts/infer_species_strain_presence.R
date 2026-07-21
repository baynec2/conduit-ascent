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

# Presence-call params (see config/snakemake.yaml). `method` selects the
# presence rule: "qvalue" (global ranked-list q-value gate) or "enrichment"
# (per-taxon per-peptide score >= margin x pooled decoy noise rate — the
# recommended default, robust at taxon scale).
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
  "Presence rule: method=%s, margin=%s, qvalue_threshold=%s, min_peptides=%s, min_confident_peptides=%s; coverage filter: score_fraction_threshold=%s, max_taxa=%s",
  method, format(margin), format(qvalue_threshold), format(min_peptides),
  format(min_confident_peptides), format(score_fraction_threshold), format(max_taxa)
)

# =============================================================================
# Read second-pass DIA-NN results
# =============================================================================
conduitR::log_with_timestamp("Reading second-pass DIA-NN parquet")

# A 0-byte parquet is the empty marker written by the second-pass DIA-NN rule
# when the second-pass database was empty (no families detected). Resolve to an
# empty detection so the run yields an empty (no-detection) conduit, rather than
# hard-failing. Guard before arrow::read_parquet, which errors on a 0-byte file.
empty_second_pass <- file.size(second_pass_diann_parquet) == 0
if (!empty_second_pass) {
  precursors <- arrow::read_parquet(second_pass_diann_parquet)
  conduitR::log_with_timestamp("Parquet rows: %d", nrow(precursors))
  empty_second_pass <- nrow(precursors) == 0
}

if (empty_second_pass) {
  conduitR::log_with_timestamp("Empty second-pass results — writing empty species/strain outputs.")
  readr::write_tsv(
    tibble::tibble(ncbi_taxonomy_id = character(), detected_taxonomy = character()),
    detected_species_strains_fp
  )
  readr::write_tsv(
    tibble::tibble(
      taxon = character(), taxon_name = character(), score = numeric(),
      n_unique_peptides_all = integer(), decoy = logical(), picked_winner = logical(),
      fdr = numeric(), qvalue = numeric(), pass = logical(),
      n_unique_peptides_q01 = integer(), score_fraction = numeric(),
      cumulative_score_fraction = numeric(), carried_forward = logical(),
      filter_reason = character()
    ),
    fdr_results_fp
  )
  sink(type = "message"); sink(); close(zz)
  quit(save = "no", status = 0)
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
# Picked target-decoy FDR at species/strain level
# =============================================================================
# conduitR::call_taxon_presence aggregates per-(taxon, decoy) scores and applies the
# picked target-decoy competition in one call (see its docs for the rationale:
# a high-abundance taxon's reversed decoy no longer outranks a true
# low-abundance taxon's target). The min_peptides gate is on
# n_unique_peptides_all.
conduitR::log_with_timestamp(
  "Running picked FDR at species/strain level (method=%s, margin=%s, qvalue_threshold=%s, min_peptides=%s)",
  method, format(margin), format(qvalue_threshold), format(min_peptides)
)

fdr_result <- conduitR::call_taxon_presence(
  pep              = psms$PEP,
  taxon            = psms$species_taxid,
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
    "WARNING: %d taxa present on only one side (target-only or decoy-only); pairing incomplete",
    fdr_result$n_missing_pair
  )
}

# Confident-peptide count per (target) taxon: number of distinct peptide
# sequences each with at least one precursor PSM at Q.Value <= 0.01. This is
# the count users typically mean when they say "N peptides per taxon"; it sits
# alongside the broader n_unique_peptides_all from call_taxon_presence (which
# includes weak-PEP contributions) so the audit trail surfaces both.
# Informational only — the presence gate (min_peptides) is on
# n_unique_peptides_all.
PSM_QVALUE_CONFIDENT <- 0.01
q01_peptide_counts <- psms |>
  dplyr::filter(!decoy, !is.na(Q.Value), Q.Value <= PSM_QVALUE_CONFIDENT) |>
  dplyr::distinct(species_taxid, Stripped.Sequence) |>
  dplyr::count(species_taxid, name = "n_unique_peptides_q01") |>
  dplyr::rename(taxon = species_taxid)

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
  "second-pass filter: %d species/strains passed picked FDR; %d carried forward (coverage mode=%s, value=%s)",
  picked_out$n_passers, picked_out$n_carried,
  picked_out$filter_mode, format(picked_out$filter_value)
)

augmented_results <- picked_out$augmented

# =============================================================================
# Format detected species/strains for output
# =============================================================================
detected_species_strains <- augmented_results |>
  dplyr::filter(carried_forward) |>
  dplyr::transmute(
    ncbi_taxonomy_id  = taxon,
    detected_taxonomy = paste0("species_strain_", taxon)
  )

conduitR::log_with_timestamp("Carried forward %d species/strains", nrow(detected_species_strains))

# =============================================================================
# Build augmented FDR results table (audit trail)
# =============================================================================
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
  # Preserve the historical column order; picked_winner + pass are the new
  # additions (the FDR computation now writes picked fdr / qvalue / pass).
  dplyr::select(
    taxon, taxon_name, score, n_unique_peptides_all, decoy,
    picked_winner, fdr, qvalue, pass, n_unique_peptides_q01,
    score_fraction, cumulative_score_fraction, carried_forward, filter_reason
  )

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
