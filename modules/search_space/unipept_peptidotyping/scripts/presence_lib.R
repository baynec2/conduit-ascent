# Pure-function helpers shared by infer_family_presence.R and
# infer_species_strain_presence.R. Kept in a separate file so the glue logic
# is reachable from testthat without the snakemake@ globals.

normalize_param <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  if (is.character(x) && (x %in% c("", "NA", "null", "None"))) return(NA_real_)
  suppressWarnings(as.numeric(x))
}

# First-pass: extract lca_taxid from Protein.Ids ("umgap|{id}|{lca_taxid}"),
# join with the taxid_map TSV to recover family_taxid for every PSM, and drop
# rows that didn't map or have a missing PEP score.
extract_family_psms <- function(precursors, taxid_map) {
  precursors |>
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
}

# Second-pass: the lca_taxid in Protein.Ids IS the species/strain taxid (the
# second-pass FASTA is built from species_strain_lca_filtered_peptides.tsv),
# so no join is needed — just extract and drop NAs.
extract_species_strain_psms <- function(precursors) {
  precursors |>
    dplyr::select(PEP, Q.Value, Stripped.Sequence, Decoy, Protein.Ids) |>
    dplyr::mutate(
      species_taxid = stringr::str_extract(Protein.Ids, "(?<=\\|)[^|]+$"),
      decoy         = as.logical(Decoy)
    ) |>
    dplyr::filter(!is.na(species_taxid), !is.na(PEP))
}

# Optional score-coverage filter + audit-table assembly for a picked-FDR pass.
# The picked target-decoy competition itself now lives in conduitR (see
# conduitR::call_taxon_presence, which aggregates PSM scores and applies the picked
# strategy in one call); this helper takes that call's per-(taxon, decoy)
# `$results` table and turns it into the schema-stable audit table the pipeline
# writes, applying the optional abundance-based coverage filter on top.
#
# Inputs:
#   picked_results — the `$results` tibble from conduitR::call_taxon_presence (cols:
#                  taxon, score, n_unique_peptides_all, decoy, picked_winner,
#                  fdr, qvalue, pass). `pass` is taken as the authoritative
#                  presence call; this helper does not re-run the FDR math.
#   q01_counts   — tibble(taxon, n_unique_peptides_q01) of confident-peptide
#                  counts (informational column carried into the audit table;
#                  NOT the gate — the min_peptides gate was on
#                  n_unique_peptides_all, inside call_taxon_presence).
#   min_peptides — the min-peptides floor used upstream; needed only to label
#                  the filter_reason of non-passing target winners (min_peptides
#                  vs fdr).
#   score_fraction_threshold / max_taxa — OPTIONAL abundance-based coverage
#                  filter applied to the picked passers (set at most one; both
#                  NA = disabled, the recall-oriented default). This is the
#                  secondary top-N gate the picked-FDR spec disables by default.
#
# Returns a list:
#   augmented      — every input row, schema-stable + score_fraction /
#                    cumulative_score_fraction / carried_forward / filter_reason.
#                    No taxon names (the caller joins those, since the lookup
#                    differs per pass).
#   detected_taxa  — character vector of carried-forward taxon ids.
#   filter_mode    — which coverage filter ran ("none (filter disabled)" etc.).
#   filter_value   — its numeric threshold (NA when disabled).
#   n_passers      — picked-FDR passers before the coverage filter.
#   n_carried      — taxa carried forward after the coverage filter.
apply_picked_presence_filter <- function(picked_results, q01_counts,
                                         min_peptides,
                                         score_fraction_threshold = NA_real_,
                                         max_taxa = NA_real_,
                                         method = "picked",
                                         min_confident_peptides = NA_real_) {
  if (!is.na(score_fraction_threshold) && !is.na(max_taxa)) {
    stop("Set at most one of score_fraction_threshold / max_taxa (both null disables the coverage filter).")
  }

  res <- picked_results |>
    dplyr::left_join(q01_counts, by = "taxon") |>
    dplyr::mutate(n_unique_peptides_q01 = ifelse(is.na(n_unique_peptides_q01),
                                                 0L, as.integer(n_unique_peptides_q01)))

  # --- Count-based presence (method = "count") -------------------------------
  # A target taxon is present iff it has at least `min_confident_peptides`
  # distinct confident (DIA-NN Q.Value <= 0.01) peptides. This bypasses the
  # picked target-decoy competition entirely: on InfiniDIA output the reported
  # decoy null is censored by --pre-filter, which corrupts any decoy-magnitude
  # statistic (enrichment) and thins the decoy-count FDR (qvalue). Confident
  # peptide COUNT, resting on DIA-NN's own (uncensored) precursor FDR, survives
  # the censoring and separates real taxa from noise. The picked fdr/qvalue
  # columns are still carried through for the audit trail but are not the gate.
  if (identical(method, "count")) {
    if (is.na(min_confident_peptides)) {
      stop("method = 'count' requires min_confident_peptides")
    }
    res <- res |>
      dplyr::mutate(
        score_fraction            = NA_real_,
        cumulative_score_fraction = NA_real_,
        carried_forward = !decoy & (n_unique_peptides_q01 >= min_confident_peptides),
        filter_reason = dplyr::case_when(
          decoy                                           ~ "decoy",
          n_unique_peptides_q01 >= min_confident_peptides ~ "",
          TRUE                                            ~ "count"
        )
      )
    detected_taxa <- res$taxon[res$carried_forward]
    return(list(
      augmented     = res,
      detected_taxa = detected_taxa,
      filter_mode   = sprintf("count (min_confident_peptides=%s)",
                              format(min_confident_peptides)),
      filter_value  = as.numeric(min_confident_peptides),
      n_passers     = length(detected_taxa),
      n_carried     = length(detected_taxa)
    ))
  }

  # Base rejection reason, before the optional coverage filter. Precedence:
  # discarded pick loser > decoy winner > passed > too-few-peptides > fails
  # q-value. Empty string means the row passed picked FDR (coverage filter
  # decides next). `pass` already encodes the q-value + min_peptides decision;
  # min_peptides is used only to split the two failure reasons.
  res <- res |>
    dplyr::mutate(
      filter_reason = dplyr::case_when(
        !picked_winner                       ~ "picked_loser",
        decoy                                ~ "decoy",
        pass                                 ~ "",
        n_unique_peptides_all < min_peptides ~ "min_peptides",
        TRUE                                 ~ "fdr"
      ),
      score_fraction            = NA_real_,
      cumulative_score_fraction = NA_real_,
      carried_forward           = FALSE
    )

  # Coverage filter applies only to the picked passers (pass == TRUE).
  pass_idx <- which(res$pass)
  pass_idx <- pass_idx[order(-res$score[pass_idx])]   # descending score
  n_passers <- length(pass_idx)

  if (n_passers == 0) {
    filter_mode  <- "none (no picked-FDR-passing taxa)"
    filter_value <- NA_real_
  } else {
    total_score <- sum(res$score[pass_idx])
    sf  <- res$score[pass_idx] / total_score
    csf <- cumsum(sf)
    res$score_fraction[pass_idx]            <- sf
    res$cumulative_score_fraction[pass_idx] <- csf

    if (!is.na(max_taxa)) {
      cutoff       <- min(as.integer(max_taxa), n_passers)
      filter_mode  <- "max_taxa"
      filter_value <- as.numeric(max_taxa)
    } else if (!is.na(score_fraction_threshold)) {
      reach        <- which(csf >= score_fraction_threshold)
      cutoff       <- if (length(reach) == 0) n_passers else reach[1]
      filter_mode  <- "score_fraction_threshold"
      filter_value <- score_fraction_threshold
    } else {
      cutoff       <- n_passers
      filter_mode  <- "none (filter disabled)"
      filter_value <- NA_real_
    }

    kept    <- pass_idx[seq_len(cutoff)]
    dropped <- if (cutoff < n_passers) pass_idx[(cutoff + 1L):n_passers] else integer(0)
    res$carried_forward[kept] <- TRUE
    if (length(dropped) > 0) res$filter_reason[dropped] <- filter_mode
  }

  detected_taxa <- res$taxon[res$carried_forward]

  list(
    augmented     = res,
    detected_taxa = detected_taxa,
    filter_mode   = filter_mode,
    filter_value  = filter_value,
    n_passers     = n_passers,
    n_carried     = length(detected_taxa)
  )
}
