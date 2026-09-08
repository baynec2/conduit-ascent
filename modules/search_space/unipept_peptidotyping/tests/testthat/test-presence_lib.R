# Integration tests for apply_picked_presence_filter() — the glue that turns a
# conduitR::call_taxon_presence picked result into the per-pass audit table. The full
# picked-FDR acceptance criteria (Poaceae rescue, 413 passers, etc.) are
# asserted in conduitR's test-call_taxon_presence.R against the real
# first_pass_fdr_results.tsv; here we check schema preservation and the
# carried_forward / filter_reason wiring.

source(testthat::test_path("..", "..", "scripts", "presence_lib.R"))

# Build a picked result the way the scripts do: call_taxon_presence on PSM-level
# input. Four taxa, each (except D) with target + decoy peptides:
#   A: strong target, weak decoy  -> target wins, passes
#   B: target wins pick but only 1 peptide -> min_peptides reject at floor 2
#   C: decoy outscores target     -> decoy wins the pick (decoy rep)
#   D: target-only, clean         -> passes
make_picked <- function(qthr = 0.5, minpep = 2) {
  psm <- tibble::tribble(
    ~taxon, ~decoy, ~peptide, ~pep,
    "A", FALSE, "A1", 0.001, "A", FALSE, "A2", 0.001, "A", FALSE, "A3", 0.002,
    "A", TRUE,  "Ad1", 0.90,
    "B", FALSE, "B1", 0.01,
    "B", TRUE,  "Bd1", 0.30, "B", TRUE, "Bd2", 0.40,
    "C", FALSE, "C1", 0.60,
    "C", TRUE,  "Cd1", 0.02, "C", TRUE, "Cd2", 0.03,
    "D", FALSE, "D1", 0.001, "D", FALSE, "D2", 0.002
  )
  conduitR::call_taxon_presence(
    pep = psm$pep, taxon = psm$taxon, decoy = psm$decoy, peptide = psm$peptide,
    qvalue_threshold = qthr, min_peptides = minpep
  )
}

q01 <- tibble::tibble(
  taxon                 = c("A", "B", "D"),
  n_unique_peptides_q01 = c(3L, 1L, 2L)
)

EXPECTED_COLS <- c(
  "taxon", "score", "n_unique_peptides_all", "decoy", "picked_winner",
  "fdr", "qvalue", "pass", "n_unique_peptides_q01",
  "score_fraction", "cumulative_score_fraction", "carried_forward",
  "filter_reason"
)

test_that("output keeps every input row and the expected columns", {
  out <- apply_picked_presence_filter(make_picked()$results, q01, min_peptides = 2)
  expect_setequal(colnames(out$augmented), EXPECTED_COLS)
  # 4 taxa -> 4 representatives; all candidate rows preserved.
  expect_equal(sum(out$augmented$picked_winner), 4L)
})

test_that("with coverage filter disabled, carried_forward == picked pass", {
  out <- apply_picked_presence_filter(make_picked()$results, q01, min_peptides = 2)
  aug <- out$augmented
  expect_equal(aug$carried_forward, aug$pass)
  passed <- aug$taxon[aug$carried_forward]
  expect_setequal(passed, c("A", "D"))
})

test_that("filter_reason vocabulary is assigned correctly", {
  out <- apply_picked_presence_filter(make_picked()$results, q01, min_peptides = 2)
  aug <- out$augmented
  reason <- function(tx, is_decoy) aug$filter_reason[aug$taxon == tx & aug$decoy == is_decoy]
  expect_equal(reason("A", FALSE), "")             # carried forward
  expect_equal(reason("A", TRUE),  "picked_loser") # A's decoy lost the pick
  expect_equal(reason("B", FALSE), "min_peptides") # target won but <2 peptides
  expect_equal(reason("C", TRUE),  "decoy")        # decoy won C's pick
  expect_equal(reason("D", FALSE), "")             # carried forward
})

test_that("min_peptides reason follows the upstream floor", {
  # With floor 1, B's single-peptide target rep passes; verify it carries.
  out <- apply_picked_presence_filter(make_picked(qthr = 0.9, minpep = 1)$results,
                                      q01, min_peptides = 1)
  expect_true("B" %in% out$augmented$taxon[out$augmented$carried_forward])
})

test_that("score-coverage max_taxa filter trims picked passers", {
  out <- apply_picked_presence_filter(make_picked()$results, q01,
                                      min_peptides = 2, max_taxa = 1)
  aug <- out$augmented
  # Only the top-scoring passer is carried forward.
  expect_equal(out$n_carried, 1L)
  # The dropped passer is tagged with the coverage filter mode.
  expect_true("max_taxa" %in% aug$filter_reason)
})

test_that("setting both coverage knobs is an error", {
  expect_error(
    apply_picked_presence_filter(make_picked()$results, q01, min_peptides = 2,
                                 score_fraction_threshold = 0.9, max_taxa = 5),
    "at most one"
  )
})
