################################################################################
# Standalone test for the logic in parse_eggnogmapper_annotations.R
#
# Tests the three failure modes found during CB023_DIA_ctrl:
#   1. Protein IDs in "tr|ACC|NAME" / "sp|ACC|NAME" format must be normalized
#   2. dplyr::na_if("-") must only apply to character columns (numeric evalue)
#   3. Column must be named "COG_category" (not the old "COG_cat")
#
# Usage: Rscript tests/scripts/test_parse_eggnogmapper_annotations.R
# Dependencies: readr, dplyr, stringr (all CRAN)
################################################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

input_file <- "tests/data/mock_emapper.emapper.annotations"

# ── Read the mock annotations (mirrors the script logic) ─────────────────────
raw <- readr::read_tsv(
  input_file,
  comment   = "##",
  col_names = TRUE,
  name_repair = "minimal",
  show_col_types = FALSE
)
colnames(raw)[1] <- "protein_id"

# ── Test 1: na_if must only apply to character columns ───────────────────────
# (was: dplyr::across(everything(), ...) — broke on numeric evalue column)
tryCatch(
  raw <- raw |>
    dplyr::mutate(dplyr::across(where(is.character), ~ dplyr::na_if(.x, "-"))),
  error = function(e) stop("FAIL test 1 (na_if on numeric): ", conditionMessage(e))
)
stopifnot("FAIL test 1: evalue must still be numeric after na_if" = is.numeric(raw$evalue))
cat("PASS test 1: na_if applied only to character columns\n")

# ── Test 2: COG_category column must exist ───────────────────────────────────
# (was: script used COG_cat — broke on real emapper 2.1.12 output)
stopifnot("FAIL test 2: COG_category column must exist" = "COG_category" %in% colnames(raw))
cat("PASS test 2: COG_category column present\n")

# ── Test 3: protein ID normalization must strip tr|/sp| prefix ───────────────
# (was: tr|P00001|PROT1_ECOLI never matched plain accessions in the join)
raw <- raw |>
  dplyr::mutate(protein_id = str_remove(protein_id, "^[a-z]+\\|") |>
                               str_remove("\\|.*$"))

stopifnot(
  "FAIL test 3: no tr|/sp| prefixes should remain" =
    !any(grepl("^(tr|sp)\\|", raw$protein_id)),
  "FAIL test 3: plain accessions P00001, P00002, A0A001 expected" =
    all(c("P00001", "P00002", "A0A001") %in% raw$protein_id)
)
cat("PASS test 3: protein IDs normalized to plain accessions\n")

cat("\nAll tests passed.\n")
