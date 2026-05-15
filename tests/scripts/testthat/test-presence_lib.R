# Tests the pure-function helpers used by infer_family_presence.R and
# infer_species_strain_presence.R. The functions live in
# modules/search_space/unipept_peptidotyping/scripts/presence_lib.R and are
# sourced from each snakemake script via snakemake@source().
#
# Scope: parsing/joining/regex glue ONLY. The FDR engine (calc_taxon_fdr) is
# tested in conduitR's own test suite; we don't retest it here.

LIB <- file.path(
  rprojroot::find_root(rprojroot::has_file("Snakefile")),
  "modules/search_space/unipept_peptidotyping/scripts/presence_lib.R"
)
source(LIB)


# ── normalize_param ──────────────────────────────────────────────────────────

testthat::test_that("normalize_param coerces nullish values to NA_real_", {
  testthat::expect_identical(normalize_param(NULL),    NA_real_)
  testthat::expect_identical(normalize_param(character(0)), NA_real_)
  testthat::expect_identical(normalize_param(""),      NA_real_)
  testthat::expect_identical(normalize_param("NA"),    NA_real_)
  testthat::expect_identical(normalize_param("null"),  NA_real_)
  testthat::expect_identical(normalize_param("None"),  NA_real_)
})

testthat::test_that("normalize_param parses numeric inputs", {
  testthat::expect_equal(normalize_param("0.9"), 0.9)
  testthat::expect_equal(normalize_param(5L),    5)
  testthat::expect_equal(normalize_param(0.01),  0.01)
})

testthat::test_that("normalize_param returns NA for unparseable strings", {
  # `suppressWarnings(as.numeric("abc"))` produces NA — important so
  # downstream `is.na()` gates trigger rather than crashing.
  testthat::expect_true(is.na(normalize_param("abc")))
})


# ── extract_family_psms ──────────────────────────────────────────────────────

make_precursors <- function(...) {
  rows <- list(...)
  do.call(rbind, lapply(rows, function(r) {
    data.frame(
      PEP               = r$PEP,
      Q.Value           = r$Q.Value           %||% 0.001,
      Stripped.Sequence = r$Stripped.Sequence %||% "AAA",
      Decoy             = r$Decoy             %||% 0L,
      Protein.Ids       = r$Protein.Ids,
      Protein.Names     = r$Protein.Names     %||% "family_E_coli",
      stringsAsFactors  = FALSE
    )
  }))
}
`%||%` <- function(a, b) if (is.null(a)) b else a

make_taxid_map <- function(...) {
  do.call(rbind, lapply(list(...), function(r) {
    data.frame(
      lca_taxid    = r$lca_taxid,
      rank         = r$rank         %||% "species",
      family_taxid = r$family_taxid,
      genus_taxid  = r$genus_taxid  %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }))
}

testthat::test_that("extract_family_psms parses umgap|id|taxid headers and joins to family", {
  precursors <- make_precursors(
    list(PEP = 0.01, Protein.Ids = "umgap|abc|562", Decoy = 0L),
    list(PEP = 0.02, Protein.Ids = "umgap|def|1280", Decoy = 1L)
  )
  taxid_map <- make_taxid_map(
    list(lca_taxid = "562",  family_taxid = "543"),
    list(lca_taxid = "1280", family_taxid = "90964")
  )
  out <- extract_family_psms(precursors, taxid_map)

  testthat::expect_equal(nrow(out), 2L)
  testthat::expect_equal(sort(out$family_taxid), c("543", "90964"))
  testthat::expect_type(out$decoy, "logical")
  testthat::expect_equal(out$decoy, c(FALSE, TRUE))
})

testthat::test_that("extract_family_psms drops PSMs whose taxid doesn't map to a family", {
  # Real-world bug class: the upstream taxid_map doesn't cover every taxid
  # observed in the DIA-NN parquet. Today those rows are silently dropped
  # by the !is.na(family_taxid) filter; this test pins that behavior.
  precursors <- make_precursors(
    list(PEP = 0.01, Protein.Ids = "umgap|abc|562"),    # mapped
    list(PEP = 0.02, Protein.Ids = "umgap|def|99999")   # NOT in map
  )
  taxid_map <- make_taxid_map(list(lca_taxid = "562", family_taxid = "543"))
  out <- extract_family_psms(precursors, taxid_map)

  testthat::expect_equal(nrow(out), 1L)
  testthat::expect_equal(out$family_taxid, "543")
})

testthat::test_that("extract_family_psms drops rows with NA PEP", {
  precursors <- make_precursors(
    list(PEP = 0.01,    Protein.Ids = "umgap|abc|562"),
    list(PEP = NA_real_, Protein.Ids = "umgap|def|562")
  )
  taxid_map <- make_taxid_map(list(lca_taxid = "562", family_taxid = "543"))
  out <- extract_family_psms(precursors, taxid_map)

  testthat::expect_equal(nrow(out), 1L)
  testthat::expect_false(any(is.na(out$PEP)))
})

testthat::test_that("extract_family_psms handles malformed Protein.Ids gracefully", {
  # A Protein.Ids string without pipes — regex (?<=\|)[^|]+$ returns NA;
  # row gets dropped via the family_taxid NA filter. No exception thrown.
  precursors <- make_precursors(
    list(PEP = 0.01, Protein.Ids = "umgap|abc|562"),   # ok
    list(PEP = 0.02, Protein.Ids = "no_pipes_here")    # malformed
  )
  taxid_map <- make_taxid_map(list(lca_taxid = "562", family_taxid = "543"))
  testthat::expect_no_error(out <- extract_family_psms(precursors, taxid_map))
  testthat::expect_equal(nrow(out), 1L)
})


# ── extract_species_strain_psms ──────────────────────────────────────────────

testthat::test_that("extract_species_strain_psms extracts the lca_taxid directly (no join)", {
  precursors <- make_precursors(
    list(PEP = 0.01, Protein.Ids = "umgap|abc|562"),
    list(PEP = 0.02, Protein.Ids = "umgap|def|511145", Decoy = 1L)
  )
  out <- extract_species_strain_psms(precursors)

  testthat::expect_equal(nrow(out), 2L)
  testthat::expect_equal(sort(out$species_taxid), c("511145", "562"))
  testthat::expect_type(out$decoy, "logical")
  testthat::expect_equal(out$decoy, c(FALSE, TRUE))
})

testthat::test_that("extract_species_strain_psms drops malformed Protein.Ids and NA PEP", {
  precursors <- make_precursors(
    list(PEP = 0.01,    Protein.Ids = "umgap|abc|562"),
    list(PEP = NA_real_, Protein.Ids = "umgap|def|562"),
    list(PEP = 0.02,    Protein.Ids = "no_pipes")
  )
  out <- extract_species_strain_psms(precursors)

  testthat::expect_equal(nrow(out), 1L)
  testthat::expect_equal(out$species_taxid, "562")
})

testthat::test_that("extract_species_strain_psms returns an empty tibble for empty input", {
  precursors <- data.frame(
    PEP               = numeric(0),
    Q.Value           = numeric(0),
    Stripped.Sequence = character(0),
    Decoy             = integer(0),
    Protein.Ids       = character(0),
    Protein.Names     = character(0),
    stringsAsFactors  = FALSE
  )
  out <- extract_species_strain_psms(precursors)
  testthat::expect_equal(nrow(out), 0L)
  testthat::expect_true(all(c("PEP", "Q.Value", "Stripped.Sequence", "decoy", "species_taxid")
                            %in% colnames(out)))
})
