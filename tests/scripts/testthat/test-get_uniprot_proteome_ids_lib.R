# Tests the organism-ID parser used by get_uniprot_proteome_ids.R. The pure
# function lives in modules/search_space/ncbi_taxonomy/scripts/get_uniprot_proteome_ids_lib.R
# and is sourced from the snakemake script via snakemake@source().

LIB <- file.path(
  rprojroot::find_root(rprojroot::has_file("Snakefile")),
  "modules/search_space/ncbi_taxonomy/scripts/get_uniprot_proteome_ids_lib.R"
)
source(LIB)


testthat::test_that("parse_organism_ids returns integers from a single-column df", {
  df <- tibble::tibble(ncbi_taxa_id = c(562L, 1280L, 1423L))
  testthat::expect_equal(parse_organism_ids(df), c(562L, 1280L, 1423L))
})

testthat::test_that("parse_organism_ids reads the first column when the df has multiple columns", {
  # Regression: the previous implementation used `readr::read_lines + as.integer`,
  # which silently coerced lines like "820\tBacteroides_uniformis_ATCC_8492" to
  # NA. The fix reads the first column directly. Real input format users have:
  df <- tibble::tibble(
    organism_id = c(820L, 562L),
    source      = c("Bacteroides_uniformis_ATCC_8492", "Escherichia_coli")
  )
  testthat::expect_equal(parse_organism_ids(df), c(820L, 562L))
})

testthat::test_that("parse_organism_ids drops non-numeric values via NA filter", {
  # If a user puts a header word or stray text in the first column, as.integer
  # returns NA; the helper filters it out rather than passing NA to UniProt.
  df <- tibble::tibble(taxa = c("562", "not_a_number", "1280"))
  testthat::expect_equal(parse_organism_ids(df), c(562L, 1280L))
})

testthat::test_that("parse_organism_ids deduplicates the result", {
  df <- tibble::tibble(taxa = c(562L, 562L, 1280L, 562L))
  testthat::expect_equal(parse_organism_ids(df), c(562L, 1280L))
})

testthat::test_that("parse_organism_ids does not append when append_id is FALSE", {
  df <- tibble::tibble(taxa = c(562L, 1280L))
  testthat::expect_equal(parse_organism_ids(df, append_id = FALSE), c(562L, 1280L))
})

testthat::test_that("parse_organism_ids appends a novel id", {
  df <- tibble::tibble(taxa = c(562L, 1280L))
  out <- parse_organism_ids(df, append_id = 9999L)
  testthat::expect_equal(out, c(562L, 1280L, 9999L))
})

testthat::test_that("parse_organism_ids does not duplicate an already-present id", {
  df <- tibble::tibble(taxa = c(562L, 1280L))
  out <- parse_organism_ids(df, append_id = 562L)
  testthat::expect_equal(out, c(562L, 1280L))
})

testthat::test_that("parse_organism_ids handles empty input cleanly", {
  # Real failure mode hit during integration: when MetaPhlAn detects no taxa
  # above threshold, the call_ncbi_taxa_ids output is header-only and the
  # parser sees a zero-row df. Returning integer(0) is the contract — the
  # crash that happens later is in conduitR::get_better_proteome_ids, not here.
  df <- tibble::tibble(taxa = integer(0))
  testthat::expect_equal(parse_organism_ids(df), integer(0))
})
