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
