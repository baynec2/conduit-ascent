# Tests the MetaPhlAn profile parser used by call_ncbi_taxa_ids.R. The pure
# function lives in modules/search_space/metaphlan/scripts/call_ncbi_taxa_ids_lib.R
# and is sourced from the snakemake script via snakemake@source().

LIB <- file.path(
  rprojroot::find_root(rprojroot::has_file("Snakefile")),
  "modules/search_space/metaphlan/scripts/call_ncbi_taxa_ids_lib.R"
)
source(LIB)

# Helper: build a MetaPhlAn-style merged-profile tibble with N samples. The
# first two columns must be exactly clade_name + NCBI_tax_id; the rest are
# numeric abundance columns named by sample.
make_profile <- function(rows, samples = "sample1") {
  df <- tibble::tibble(
    clade_name  = vapply(rows, `[[`, character(1), "clade_name"),
    NCBI_tax_id = vapply(rows, `[[`, character(1), "NCBI_tax_id")
  )
  for (s in samples) {
    df[[s]] <- vapply(rows, function(r) {
      if (is.null(r[[s]])) 0 else as.numeric(r[[s]])
    }, numeric(1))
  }
  df
}

# Canonical MetaPhlAn rows used as building blocks.
ECOLI <- list(
  clade_name  = "k__Bacteria|p__Pseudomonadota|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Escherichia|s__Escherichia_coli",
  NCBI_tax_id = "2|1224|1236|91347|543|561|562"
)
BSUB <- list(
  clade_name  = "k__Bacteria|p__Bacillota|c__Bacilli|o__Bacillales|f__Bacillaceae|g__Bacillus|s__Bacillus_subtilis",
  NCBI_tax_id = "2|1239|91061|1385|186817|1386|1423"
)
# Partial-taxonomy row — only resolved to family. species_ncbi should be NA
# after the separate(fill = "right") and get filtered out.
FAMILY_ONLY <- list(
  clade_name  = "k__Bacteria|p__Pseudomonadota|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae",
  NCBI_tax_id = "2|1224|1236|91347|543"
)


testthat::test_that("parse_metaphlan_profiles returns species_ncbi for taxa above threshold", {
  profiles <- make_profile(list(c(ECOLI, sample1 = 50.0)))
  out <- parse_metaphlan_profiles(profiles, threshold = 0.01)

  testthat::expect_s3_class(out, "tbl_df")
  testthat::expect_named(out, "ncbi_taxonomy_id")
  testthat::expect_equal(out$ncbi_taxonomy_id, "562")
})

testthat::test_that("parse_metaphlan_profiles drops taxa below threshold", {
  profiles <- make_profile(list(
    c(ECOLI, sample1 = 50.0),
    c(BSUB,  sample1 = 0.005)   # below 0.01 threshold
  ))
  out <- parse_metaphlan_profiles(profiles, threshold = 0.01)

  testthat::expect_equal(out$ncbi_taxonomy_id, "562")
})

testthat::test_that("parse_metaphlan_profiles deduplicates across samples", {
  # Same species above threshold in both samples → should appear once.
  profiles <- make_profile(
    list(c(ECOLI, sample1 = 50.0, sample2 = 30.0)),
    samples = c("sample1", "sample2")
  )
  out <- parse_metaphlan_profiles(profiles, threshold = 0.01)

  testthat::expect_equal(nrow(out), 1L)
  testthat::expect_equal(out$ncbi_taxonomy_id, "562")
})

testthat::test_that("parse_metaphlan_profiles keeps species above threshold in any sample", {
  # E. coli above threshold in sample1 only; B. subtilis above only in sample2.
  # Both should appear in the result (OR semantics across samples).
  profiles <- make_profile(
    list(
      c(ECOLI, sample1 = 50.0, sample2 = 0.001),
      c(BSUB,  sample1 = 0.001, sample2 = 30.0)
    ),
    samples = c("sample1", "sample2")
  )
  out <- parse_metaphlan_profiles(profiles, threshold = 0.01)

  testthat::expect_setequal(out$ncbi_taxonomy_id, c("562", "1423"))
})

testthat::test_that("parse_metaphlan_profiles drops rows without species_ncbi", {
  # A row resolved only to family (no s__ rank) yields NA species_ncbi after
  # the separate(fill = "right") and must be excluded.
  profiles <- make_profile(list(
    c(ECOLI,        sample1 = 50.0),
    c(FAMILY_ONLY,  sample1 = 50.0)
  ))
  out <- parse_metaphlan_profiles(profiles, threshold = 0.01)

  testthat::expect_equal(out$ncbi_taxonomy_id, "562")
})

testthat::test_that("parse_metaphlan_profiles returns empty result when all below threshold", {
  # Real-world failure mode: 1k-read FASTQ produced 100% UNCLASSIFIED, so
  # call_ncbi_taxa_ids emits only a header. The pure function must return a
  # zero-row tibble that still has the ncbi_taxonomy_id column so downstream
  # readers don't choke on the schema.
  profiles <- make_profile(list(c(ECOLI, sample1 = 0.001)))
  out <- parse_metaphlan_profiles(profiles, threshold = 0.01)

  testthat::expect_equal(nrow(out), 0L)
  testthat::expect_named(out, "ncbi_taxonomy_id")
})
